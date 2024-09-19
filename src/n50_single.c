/*
 * N50 Calculator v2.0
 * Andrea Telatin, 2024
 *
 * This program calculates N50 and other sequence statistics from FASTA or FASTQ files.
 *
 * Features:
 * - Supports FASTA and FASTQ formats (including gzipped files)
 * - Multi-threaded processing for improved performance
 * - Automatic format detection based on file extension
 * - Optional header output and N50-only output
 *
 * Usage: n50 [options] [filename]
 * Options:
 *   --fasta/-a: Force FASTA input
 *   --fastq/-q: Force FASTQ input
 *   --header/-h: Print header in output
 *   --n50/-n: Output only N50 value
 *
 * Input: File or stdin
 * Output: Tab-separated values (Format, Total Length, Total Sequences, N50)
 *
 * Compile: gcc -o n50 n50.c -lz -lpthread -O3
 *
 * Dependencies: zlib, pthread
 *
 * Note: Max threads set to 8, initial capacity 1,000,000 sequences.
 *       Adjust MAX_THREADS and INITIAL_CAPACITY for different requirements.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <pthread.h>
#include <ctype.h>
#include <getopt.h>
#include <stdbool.h>

#define BUFFER_SIZE 1024 * 1024  // 1MB buffer
#define MAX_THREADS 8
#define INITIAL_CAPACITY 1000000

typedef struct {
    long long *lengths;
    int start;
    int end;
} ThreadData;

long long total_length = 0;
int *lengths;
int length_count = 0;
int length_capacity = INITIAL_CAPACITY;
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

int compare(const void *a, const void *b) {
    return (*(int*)b - *(int*)a);
}

void *process_chunk(void *arg) {
    ThreadData *data = (ThreadData*)arg;
    for (int i = data->start; i < data->end; i++) {
        pthread_mutex_lock(&mutex);
        total_length += data->lengths[i];
        if (length_count == length_capacity) {
            length_capacity *= 2;
            lengths = realloc(lengths, length_capacity * sizeof(int));
        }
        lengths[length_count++] = data->lengths[i];
        pthread_mutex_unlock(&mutex);
    }
    return NULL;
}

void process_fasta(gzFile fp) {
    char buffer[BUFFER_SIZE];
    long long *chunk_lengths = NULL;
    int chunk_count = 0;
    int chunk_capacity = INITIAL_CAPACITY;
    int current_length = 0;

    chunk_lengths = malloc(chunk_capacity * sizeof(long long));

    while (gzgets(fp, buffer, BUFFER_SIZE) != NULL) {
        if (buffer[0] == '>') {
            if (current_length > 0) {
                if (chunk_count == chunk_capacity) {
                    chunk_capacity *= 2;
                    chunk_lengths = realloc(chunk_lengths, chunk_capacity * sizeof(long long));
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
            chunk_lengths = realloc(chunk_lengths, chunk_capacity * sizeof(long long));
        }
        chunk_lengths[chunk_count++] = current_length;
    }

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

    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    free(chunk_lengths);
}

void process_fastq(gzFile fp) {
    char buffer[BUFFER_SIZE];
    int line_count = 0;
    long long *chunk_lengths = NULL;
    int chunk_count = 0;
    int chunk_capacity = INITIAL_CAPACITY;

    chunk_lengths = malloc(chunk_capacity * sizeof(long long));

    while (gzgets(fp, buffer, BUFFER_SIZE) != NULL) {
        line_count++;
        if (line_count % 4 == 2) {  // Sequence line
            int len = strcspn(buffer, "\n");  // Count characters until newline
            if (chunk_count == chunk_capacity) {
                chunk_capacity *= 2;
                chunk_lengths = realloc(chunk_lengths, chunk_capacity * sizeof(long long));
            }
            chunk_lengths[chunk_count++] = len;
        }
    }

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

    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    free(chunk_lengths);
}

bool is_fastq_filename(const char *filename) {
    const char *ext = strrchr(filename, '.');
    if (ext != NULL) {
        if (strcmp(ext, ".fastq") == 0 || strcmp(ext, ".fq") == 0) {
            return true;
        }
        if (strcmp(ext, ".gz") == 0) {
            const char *prev_ext = ext - 1;
            while (prev_ext > filename && *prev_ext != '.') {
                prev_ext--;
            }
            if (strcmp(prev_ext, ".fastq.gz") == 0 || strcmp(prev_ext, ".fq.gz") == 0) {
                return true;
            }
        }
    }
    return false;
}

int main(int argc, char *argv[]) {
      gzFile fp;
    bool is_fastq = false;
    bool opt_header = false;
    bool opt_n50 = false;
    char *filename = NULL;
    bool format_specified = false;

    static struct option long_options[] = {
        {"fasta", no_argument, 0, 'a'},
        {"fastq", no_argument, 0, 'q'},
        {"header", no_argument, 0, 'h'},
        {"n50", no_argument, 0, 'n'},
        {0, 0, 0, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "aqhn", long_options, NULL)) != -1) {
        switch (opt) {
            case 'a':
                is_fastq = false;
                format_specified = true;
                break;
            case 'q':
                is_fastq = true;
                format_specified = true;
                break;
            case 'h':
                opt_header = true;
                break;
            case 'n':
                opt_n50 = true;
                break;
            default:
                fprintf(stderr, "Usage: %s [--fasta | --fastq] [--header] [--n50] [filename]\n", argv[0]);
                return 1;
        }
    }

    if (optind < argc) {
        filename = argv[optind];
        if (!format_specified) {
            is_fastq = is_fastq_filename(filename);
        }
    }

    if (filename) {
        fp = gzopen(filename, "r");
        if (!fp) {
            fprintf(stderr, "Error: Cannot open file %s\n", filename);
            return 1;
        }
    } else {
        fp = gzdopen(fileno(stdin), "r");
    }

    lengths = malloc(length_capacity * sizeof(int));

    if (is_fastq) {
        process_fastq(fp);
    } else {
        process_fasta(fp);
    }

    gzclose(fp);

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
    
    if (opt_header && !opt_n50) {
        printf("Format\tTotal_Length\tTotal_Sequences\tN50\n");
    }
    if (opt_n50) {
        printf("%d\n", n50);
    } else {
        printf("%s\t%lld\t%d\t%d\n", is_fastq ? "FASTQ" : "FASTA", total_length, length_count, n50);
    }

    free(lengths);
    return 0;
}
