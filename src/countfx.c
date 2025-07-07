#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <zlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>

#define CHUNK_SIZE 1048576 // 1MB chunks
#define NUM_THREADS 4

typedef enum {
    FORMAT_AUTO,
    FORMAT_FASTA,
    FORMAT_FASTQ
} FileFormat;

typedef struct {
    gzFile file;
    size_t count;
    int thread_id;
    FileFormat format;
} ThreadData;

pthread_mutex_t file_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t count_mutex = PTHREAD_MUTEX_INITIALIZER;
int global_error = 0;

FileFormat detect_format_from_extension(const char *filename) {
    const char *ext = strrchr(filename, '.');
    if (!ext) return FORMAT_AUTO;
    
    if (strcmp(ext, ".gz") == 0) {
        char *temp = strdup(filename);
        temp[strlen(temp) - 3] = '\0';
        ext = strrchr(temp, '.');
        if (ext) {
            if (strcmp(ext, ".fq") == 0 || strcmp(ext, ".fastq") == 0) {
                free(temp);
                return FORMAT_FASTQ;
            }
        }
        free(temp);
        return FORMAT_FASTA;
    }
    
    if (strcmp(ext, ".fq") == 0 || strcmp(ext, ".fastq") == 0) {
        return FORMAT_FASTQ;
    }
    
    return FORMAT_FASTA;
}

void *count_sequences(void *arg) {
    ThreadData *data = (ThreadData *)arg;
    char *buffer = NULL;
    int bytes_read;
    size_t local_count = 0;

    buffer = (char *)aligned_alloc(8, CHUNK_SIZE);
    if (buffer == NULL) {
        fprintf(stderr, "Thread %d: Memory allocation failed\n", data->thread_id);
        global_error = 1;
        return NULL;
    }

    while (!global_error) {
        pthread_mutex_lock(&file_mutex);
        bytes_read = gzread(data->file, buffer, CHUNK_SIZE);
        int gz_err;
        const char *gz_err_msg = gzerror(data->file, &gz_err);
        pthread_mutex_unlock(&file_mutex);

        if (bytes_read < 0 || gz_err != Z_OK) {
            fprintf(stderr, "Thread %d: gzread error: %s\n", data->thread_id, gz_err_msg);
            global_error = 1;
            break;
        }

        if (bytes_read == 0) break; // End of file

        if (data->format == FORMAT_FASTA) {
            for (int i = 0; i < bytes_read; i++) {
                if (buffer[i] == '>') {
                    local_count++;
                }
            }
        } else if (data->format == FORMAT_FASTQ) {
            for (int i = 0; i < bytes_read; i++) {
                if (buffer[i] == '\n') {
                    local_count++;
                }
            }
        }
    }

    pthread_mutex_lock(&count_mutex);
    data->count += local_count;
    pthread_mutex_unlock(&count_mutex);

    free(buffer);
    return NULL;
}

void print_usage(const char *program_name) {
    fprintf(stderr, "Usage: %s [options] <file>\n", program_name);
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  --fasta    Force FASTA format\n");
    fprintf(stderr, "  --fastq    Force FASTQ format\n");
    fprintf(stderr, "  -h, --help Show this help message\n");
    fprintf(stderr, "\nFile format detection:\n");
    fprintf(stderr, "  .fq, .fq.gz, .fastq, .fastq.gz -> FASTQ\n");
    fprintf(stderr, "  All other extensions -> FASTA\n");
}

int main(int argc, char **argv) {
    FileFormat format = FORMAT_AUTO;
    char *filename = NULL;
    
    static struct option long_options[] = {
        {"fasta", no_argument, 0, 'a'},
        {"fastq", no_argument, 0, 'q'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };
    
    int option_index = 0;
    int c;
    
    while ((c = getopt_long(argc, argv, "aqh", long_options, &option_index)) != -1) {
        switch (c) {
            case 'a':
                format = FORMAT_FASTA;
                break;
            case 'q':
                format = FORMAT_FASTQ;
                break;
            case 'h':
                print_usage(argv[0]);
                return 0;
            case '?':
                print_usage(argv[0]);
                return 1;
            default:
                print_usage(argv[0]);
                return 1;
        }
    }
    
    if (optind >= argc) {
        fprintf(stderr, "Error: No input file specified\n");
        print_usage(argv[0]);
        return 1;
    }
    
    filename = argv[optind];
    
    if (format == FORMAT_AUTO) {
        format = detect_format_from_extension(filename);
    }

    gzFile file = gzopen(filename, "rb");
    if (!file) {
        fprintf(stderr, "Error: Cannot open file %s: %s\n", filename, strerror(errno));
        return 1;
    }

    pthread_t threads[NUM_THREADS];
    ThreadData thread_data[NUM_THREADS];

    for (int i = 0; i < NUM_THREADS; i++) {
        thread_data[i].file = file;
        thread_data[i].count = 0;
        thread_data[i].thread_id = i;
        thread_data[i].format = format;

        if (pthread_create(&threads[i], NULL, count_sequences, &thread_data[i]) != 0) {
            fprintf(stderr, "Error creating thread %d: %s\n", i, strerror(errno));
            gzclose(file);
            return 1;
        }
    }

    for (int i = 0; i < NUM_THREADS; i++) {
        pthread_join(threads[i], NULL);
    }

    gzclose(file);

    if (global_error) {
        fprintf(stderr, "An error occurred during processing. Results may be incomplete.\n");
        return 1;
    }

    size_t total_count = 0;
    for (int i = 0; i < NUM_THREADS; i++) {
        total_count += thread_data[i].count;
    }

    size_t sequence_count;
    if (format == FORMAT_FASTQ) {
        sequence_count = total_count / 4;
    } else {
        sequence_count = total_count;
    }

    printf("Total sequences: %zu\n", sequence_count);
    return 0;
}