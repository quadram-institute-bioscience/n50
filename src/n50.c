#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <pthread.h>
#include <ctype.h>
#include <getopt.h>
#include <stdbool.h>

#define BUFFER_SIZE 4 * 1024 * 1024  // 1MB buffer
#define MAX_THREADS 8
#define INITIAL_CAPACITY 1000000
#define VERSION "1.0.0"

typedef struct {
    int *lengths;
    int start;
    int end;
    long long total_length;
} ThreadData;

typedef struct {
    char *filename;
    long long total_length;
    int length_count;
    int n50;
    bool is_fastq;
} FileStats;

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

void *process_chunk(void *arg) {
    ThreadData *data = (ThreadData*)arg;
    long long local_total_length = 0;
    int *local_lengths = malloc((data->end - data->start) * sizeof(int));
    int local_count = 0;

    for (int i = data->start; i < data->end; i++) {
        local_total_length += data->lengths[i];
        local_lengths[local_count++] = data->lengths[i];
    }

    ThreadData *result = malloc(sizeof(ThreadData));
    result->lengths = local_lengths;
    result->start = data->start;
    result->end = data->end;
    result->total_length = local_total_length;

    return result;
}

FileStats process_file(const char *filename, bool force_fasta, bool force_fastq) {
    gzFile fp = gzopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open file %s\n", filename);
        exit(1);
    }

    char buffer[BUFFER_SIZE];
    int *chunk_lengths = NULL;
    int chunk_count = 0;
    int chunk_capacity = INITIAL_CAPACITY;
    int current_length = 0;
    bool is_fastq = force_fastq || (!force_fasta && (strstr(filename, ".fastq") != NULL || strstr(filename, ".fq") != NULL));

    chunk_lengths = malloc(chunk_capacity * sizeof(int));

    if (is_fastq) {
        int line_count = 0;
        while (gzgets(fp, buffer, BUFFER_SIZE) != NULL) {
            line_count++;
            if (line_count % 4 == 2) {  // Sequence line
                int len = strcspn(buffer, "\n");  // Count characters until newline
                if (chunk_count == chunk_capacity) {
                    chunk_capacity *= 2;
                    chunk_lengths = realloc(chunk_lengths, chunk_capacity * sizeof(int));
                }
                chunk_lengths[chunk_count++] = len;
            }
        }
    } else {
        while (gzgets(fp, buffer, BUFFER_SIZE) != NULL) {
            if (buffer[0] == '>') {
                if (current_length > 0) {
                    if (chunk_count == chunk_capacity) {
                        chunk_capacity *= 2;
                        chunk_lengths = realloc(chunk_lengths, chunk_capacity * sizeof(int));
                    }
                    chunk_lengths[chunk_count++] = current_length;
                    current_length = 0;
                }
            } else {
                for (char *p = buffer; *p != '\0'; p++) {
                    if (isalpha(*p)) {
                        current_length++;
                    }
                }
            }
        }
        if (current_length > 0) {
            if (chunk_count == chunk_capacity) {
                chunk_capacity *= 2;
                chunk_lengths = realloc(chunk_lengths, chunk_capacity * sizeof(int));
            }
            chunk_lengths[chunk_count++] = current_length;
        }
    }

    gzclose(fp);

    // Process chunks using threads
    int num_threads = (chunk_count < MAX_THREADS) ? chunk_count : MAX_THREADS;
    pthread_t threads[MAX_THREADS];
    ThreadData thread_data[MAX_THREADS];
    int chunk_size = chunk_count / num_threads;

    for (int i = 0; i < num_threads; i++) {
        thread_data[i].lengths = chunk_lengths;
        thread_data[i].start = i * chunk_size;
        thread_data[i].end = (i == num_threads - 1) ? chunk_count : (i + 1) * chunk_size;
        pthread_create(&threads[i], NULL, process_chunk, &thread_data[i]);
    }

    long long total_length = 0;
    int *lengths = malloc(chunk_count * sizeof(int));
    int length_count = 0;

    for (int i = 0; i < num_threads; i++) {
        ThreadData *result;
        pthread_join(threads[i], (void**)&result);
        total_length += result->total_length;
        memcpy(lengths + length_count, result->lengths, (result->end - result->start) * sizeof(int));
        length_count += result->end - result->start;
        free(result->lengths);
        free(result);
    }

    free(chunk_lengths);

    // Sort lengths in descending order
    qsort(lengths, length_count, sizeof(int), compare);

    // Calculate N50
    long long cumulative_length = 0;
    long long half_total_length = total_length / 2;
    int n50 = 0;
    for (int i = 0; i < length_count; i++) {
        cumulative_length += lengths[i];
        if (cumulative_length >= half_total_length) {
            n50 = lengths[i];
            break;
        }
    }

    free(lengths);

    FileStats stats = {
        .filename = strdup(filename),
        .total_length = total_length,
        .length_count = length_count,
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
