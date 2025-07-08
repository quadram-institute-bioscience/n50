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
#define DEFAULT_THREADS 4

static int num_threads = DEFAULT_THREADS;
typedef struct {
    gzFile file;
    size_t seq_count;
    int thread_id;
} ThreadData;
pthread_mutex_t file_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t count_mutex = PTHREAD_MUTEX_INITIALIZER;
int global_error = 0;
void *count_newlines(void *arg) {
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
        for (int i = 0; i < bytes_read; i++) {
            if (buffer[i] == '>') {
                local_count++;
            }
        }
    }
    pthread_mutex_lock(&count_mutex);
    data->seq_count += local_count;
    pthread_mutex_unlock(&count_mutex);
    free(buffer);
    return NULL;
}
int main(int argc, char **argv) {
    int opt;
    static struct option long_options[] = {
        {"threads", required_argument, 0, 't'},
        {0, 0, 0, 0}
    };

    while ((opt = getopt_long(argc, argv, "t:", long_options, NULL)) != -1) {
        switch (opt) {
            case 't':
                num_threads = atoi(optarg);
                if (num_threads <= 0) {
                    fprintf(stderr, "Invalid thread count: %s\n", optarg);
                    return 1;
                }
                break;
            default:
                fprintf(stderr, "Usage: %s [--threads N] <fasta.gz file>\n", argv[0]);
                return 1;
        }
    }

    if (optind >= argc) {
        fprintf(stderr, "Usage: %s [--threads N] <fasta.gz file>\n", argv[0]);
        return 1;
    }

    const char *filename = argv[optind];
    gzFile file;
    if (strcmp(filename, "-") == 0) {
        file = gzdopen(fileno(stdin), "rb");
    } else {
        file = gzopen(filename, "rb");
    }
    if (!file) {
        fprintf(stderr, "Error: Cannot open file %s: %s\n", filename, strerror(errno));
        return 1;
    }

    int first = gzgetc(file);
    if (first == -1) {
        fprintf(stderr, "Error reading from %s\n", filename);
        gzclose(file);
        return 1;
    }
    if (first != '>') {
        fprintf(stderr, "Invalid FASTA format: first character is not '>'\n");
        gzclose(file);
        return 1;
    }
    gzungetc(first, file);

    pthread_t *threads = malloc(sizeof(pthread_t) * num_threads);
    ThreadData *thread_data = malloc(sizeof(ThreadData) * num_threads);
    if (!threads || !thread_data) {
        fprintf(stderr, "Memory allocation failure\n");
        gzclose(file);
        free(threads);
        free(thread_data);
        return 1;
    }
    for (int i = 0; i < num_threads; i++) {
        thread_data[i].file = file;
        thread_data[i].seq_count = 0;
        thread_data[i].thread_id = i;
        if (pthread_create(&threads[i], NULL, count_newlines, &thread_data[i]) != 0) {
            fprintf(stderr, "Error creating thread %d: %s\n", i, strerror(errno));
            gzclose(file);
            free(threads);
            free(thread_data);
            return 1;
        }
    }
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }
    gzclose(file);
    if (global_error) {
        fprintf(stderr, "An error occurred during processing. Results may be incomplete.\n");
        free(threads);
        free(thread_data);
        return 1;
    }
    size_t total_sequences = 0;
    for (int i = 0; i < num_threads; i++) {
        total_sequences += thread_data[i].newline_count;
    }
    printf("Total sequences: %zu\n", total_newlines);
    free(threads);
    free(thread_data);
    return 0;
}