#define GNUSOURCE
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <zlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#define CHUNK_SIZE 1048576 // 1MB chunks
#define NUM_THREADS 4
typedef struct {
    gzFile file;
    size_t newline_count;
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
    data->newline_count += local_count;
    pthread_mutex_unlock(&count_mutex);
    free(buffer);
    return NULL;
}
int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <fasta.gz file>\n", argv[0]);
        return 1;
    }
    gzFile file = gzopen(argv[1], "rb");
    if (!file) {
        fprintf(stderr, "Error: Cannot open file %s: %s\n", argv[1], strerror(errno));
        return 1;
    }
    pthread_t threads[NUM_THREADS];
    ThreadData thread_data[NUM_THREADS];
    for (int i = 0; i < NUM_THREADS; i++) {
        thread_data[i].file = file;
        thread_data[i].newline_count = 0;
        thread_data[i].thread_id = i;
        if (pthread_create(&threads[i], NULL, count_newlines, &thread_data[i]) != 0) {
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
    size_t total_newlines = 0;
    for (int i = 0; i < NUM_THREADS; i++) {
        total_newlines += thread_data[i].newline_count;
    }
    printf("Total sequences: %zu\n", total_newlines);
    return 0;
}