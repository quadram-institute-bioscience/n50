#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <ctype.h>
#include <getopt.h>
#include <stdbool.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>

#define MAX_THREADS 16
#define INITIAL_CAPACITY 1000000
#define VERSION "2.0.0"

typedef struct {
    char *filename;
    long long total_length;
    int length_count;
    int n50;
    bool is_fastq;
} FileStats;

typedef struct {
    char *start;
    size_t size;
    int *lengths;
    int capacity;
    int count;
    long long total_length;
    bool is_fastq;
} ThreadData;

void print_usage(const char *program_name) {
    printf("Usage: %s [OPTIONS] [FILENAME...]\n", program_name);
    printf("\nOptions:\n");
    printf("  -a, --fasta        Force FASTA input format\n");
    printf("  -q, --fastq        Force FASTQ input format\n");
    printf("  -h, --header       Print header in output\n");
    printf("  -n, --n50          Output only N50 value\n");
    printf("      --help         Display this help message and exit\n");
    printf("      --version      Display version information and exit\n");
    printf("\nDescription:\n");
    printf("  Calculate N50 and other sequence statistics from FASTA or FASTQ files.\n");
    printf("  Supports multiple input files and automatic format detection.\n");
}

void print_version() {
    printf("N50 Calculator version %s\n", VERSION);
    printf("Copyright (C) 2024 Andrea Telatin\n");
    printf("License: MIT\n");
}

int compare(const void *a, const void *b) {
    return (*(int*)b - *(int*)a);
}

void *parallel_read_and_parse(void *arg) {
    ThreadData *data = (ThreadData *)arg;
    char *current = data->start;
    char *end = data->start + data->size;
    int current_length = 0;

    while (current < end) {
        if (data->is_fastq) {
            // FASTQ parsing
            if (*current == '@') {
                while (current < end && *current != '\n') current++;
                current++;  // Skip newline
                char *seq_start = current;
                while (current < end && *current != '\n') current++;
                int seq_length = current - seq_start;
                
                if (data->count == data->capacity) {
                    data->capacity *= 2;
                    data->lengths = realloc(data->lengths, data->capacity * sizeof(int));
                }
                data->lengths[data->count++] = seq_length;
                data->total_length += seq_length;

                // Skip quality line
                current += 2;  // Skip '+' and newline
                while (current < end && *current != '\n') current++;
                current++;  // Skip newline
            } else {
                current++;
            }
        } else {
            // FASTA parsing
            if (*current == '>') {
                if (current_length > 0) {
                    if (data->count == data->capacity) {
                        data->capacity *= 2;
                        data->lengths = realloc(data->lengths, data->capacity * sizeof(int));
                    }
                    data->lengths[data->count++] = current_length;
                    data->total_length += current_length;
                    current_length = 0;
                }
                while (current < end && *current != '\n') current++;
                current++;  // Skip newline
            } else {
                current_length += (*current >= 'A' && *current <= 'Z') || (*current >= 'a' && *current <= 'z');
                current++;
            }
        }
    }

    // Handle last sequence for FASTA
    if (!data->is_fastq && current_length > 0) {
        if (data->count == data->capacity) {
            data->capacity *= 2;
            data->lengths = realloc(data->lengths, data->capacity * sizeof(int));
        }
        data->lengths[data->count++] = current_length;
        data->total_length += current_length;
    }

    return NULL;
}

void parallel_merge_sort(int *arr, int *temp, int left, int right) {
    if (left < right) {
        int mid = left + (right - left) / 2;
        #pragma omp task
        parallel_merge_sort(arr, temp, left, mid);
        #pragma omp task
        parallel_merge_sort(arr, temp, mid + 1, right);
        #pragma omp taskwait

        // Merge
        int i = left, j = mid + 1, k = left;
        while (i <= mid && j <= right) {
            if (arr[i] >= arr[j]) {
                temp[k++] = arr[i++];
            } else {
                temp[k++] = arr[j++];
            }
        }
        while (i <= mid) temp[k++] = arr[i++];
        while (j <= right) temp[k++] = arr[j++];
        for (i = left; i <= right; i++) arr[i] = temp[i];
    }
}

FileStats process_file(const char *filename, bool force_fasta, bool force_fastq) {
    int fd = open(filename, O_RDONLY);
    if (fd == -1) {
        fprintf(stderr, "Error: Cannot open file %s\n", filename);
        exit(1);
    }

    off_t file_size = lseek(fd, 0, SEEK_END);
    char *file_content = mmap(NULL, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
    if (file_content == MAP_FAILED) {
        fprintf(stderr, "Error: Cannot mmap file %s\n", filename);
        close(fd);
        exit(1);
    }

    bool is_fastq = force_fastq || (!force_fasta && (strstr(filename, ".fastq") != NULL || strstr(filename, ".fq") != NULL));
    int num_threads = MAX_THREADS;
    size_t chunk_size = file_size / num_threads;

    ThreadData thread_data[MAX_THREADS];
    pthread_t threads[MAX_THREADS];

    for (int i = 0; i < num_threads; i++) {
        thread_data[i].start = file_content + i * chunk_size;
        thread_data[i].size = (i == num_threads - 1) ? file_size - i * chunk_size : chunk_size;
        thread_data[i].lengths = malloc(INITIAL_CAPACITY * sizeof(int));
        thread_data[i].capacity = INITIAL_CAPACITY;
        thread_data[i].count = 0;
        thread_data[i].total_length = 0;
        thread_data[i].is_fastq = is_fastq;
        pthread_create(&threads[i], NULL, parallel_read_and_parse, &thread_data[i]);
    }

    long long total_length = 0;
    int total_count = 0;
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
        total_length += thread_data[i].total_length;
        total_count += thread_data[i].count;
    }

    // Merge sorted arrays
    int *merged_lengths = malloc(total_count * sizeof(int));
    int *temp = malloc(total_count * sizeof(int));
    int index = 0;
    for (int i = 0; i < num_threads; i++) {
        memcpy(merged_lengths + index, thread_data[i].lengths, thread_data[i].count * sizeof(int));
        index += thread_data[i].count;
        free(thread_data[i].lengths);
    }

    // Parallel sort
    #pragma omp parallel
    {
        #pragma omp single
        parallel_merge_sort(merged_lengths, temp, 0, total_count - 1);
    }

    free(temp);

    // Calculate N50
    long long cumulative_length = 0;
    long long half_total_length = total_length / 2;
    int n50 = 0;
    #pragma omp parallel for reduction(+:cumulative_length) shared(n50)
    for (int i = 0; i < total_count; i++) {
        cumulative_length += merged_lengths[i];
        if (cumulative_length >= half_total_length && n50 == 0) {
            #pragma omp critical
            {
                if (n50 == 0) n50 = merged_lengths[i];
            }
        }
    }

    free(merged_lengths);
    munmap(file_content, file_size);
    close(fd);

    FileStats stats = {
        .filename = strdup(filename),
        .total_length = total_length,
        .length_count = total_count,
        .n50 = n50,
        .is_fastq = is_fastq
    };

    return stats;
}

int main(int argc, char *argv[]) {
    bool opt_header = false;
    bool opt_n50 = false;
    bool force_fasta = false;
    bool force_fastq = false;

    static struct option long_options[] = {
        {"fasta", no_argument, 0, 'a'},
        {"fastq", no_argument, 0, 'q'},
        {"header", no_argument, 0, 'H'},
        {"n50", no_argument, 0, 'n'},
        {"help", no_argument, 0, 'h'},
        {"version", no_argument, 0, 'v'},
        {0, 0, 0, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "aqhnHv", long_options, NULL)) != -1) {
        switch (opt) {
            case 'a':
                force_fasta = true;
                break;
            case 'q':
                force_fastq = true;
                break;
            case 'H':
                opt_header = true;
                break;
            case 'n':
                opt_n50 = true;
                break;
            case 'h':
                print_usage(argv[0]);
                return 0;
            case 'v':
                print_version();
                return 0;
            default:
                fprintf(stderr, "Try '%s --help' for more information.\n", argv[0]);
                return 1;
        }
    }

    if (optind == argc) {
        fprintf(stderr, "Error: No input files specified\n");
        fprintf(stderr, "Try '%s --help' for more information.\n", argv[0]);
        return 1;
    }

    if (opt_header && !opt_n50) {
        printf("Filename\tFormat\tTotal_Length\tTotal_Sequences\tN50\n");
    }

    for (int i = optind; i < argc; i++) {
        FileStats stats = process_file(argv[i], force_fasta, force_fastq);
        
        if (opt_n50) {
            printf("%s\t%d\n", stats.filename, stats.n50);
        } else {
            printf("%s\t%s\t%lld\t%d\t%d\n", stats.filename, stats.is_fastq ? "FASTQ" : "FASTA", stats.total_length, stats.length_count, stats.n50);
        }

        free(stats.filename);
    }

    return 0;
}