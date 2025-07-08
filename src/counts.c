#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <zlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <stdatomic.h>

#define CHUNK_SIZE 1048576 // 1MB chunks
#define NUM_THREADS 4
#define MAX_THREADS 64

typedef struct {
    char *buffer;
    int bytes_read;
} Chunk;

typedef struct ChunkNode {
    Chunk chunk;
    struct ChunkNode *next;
} ChunkNode;

typedef struct {
    ChunkNode *head;
    ChunkNode *tail;
    pthread_mutex_t mutex;
    pthread_cond_t cond;
    int producer_done;
} ChunkQueue;

typedef struct {
    ChunkQueue *queue;
    size_t newline_count;
    int thread_id;
} ThreadData;

pthread_mutex_t file_mutex = PTHREAD_MUTEX_INITIALIZER;
atomic_int global_error = 0;

void *read_chunks(void *arg) {
    gzFile file = ((void **)arg)[0];
    ChunkQueue *queue = ((void **)arg)[1];
    char *buffer;
    int bytes_read;
    int gz_err;
    const char *gz_err_msg;

    while (1) {
        buffer = (char *)malloc(CHUNK_SIZE);
        if (buffer == NULL) {
            fprintf(stderr, "Producer: Memory allocation failed\n");
            atomic_store(&global_error, 1);
            break;
        }

        pthread_mutex_lock(&file_mutex);
        bytes_read = gzread(file, buffer, CHUNK_SIZE);
        gz_err_msg = gzerror(file, &gz_err);
        pthread_mutex_unlock(&file_mutex);

        if (bytes_read < 0 || gz_err != Z_OK) {
            fprintf(stderr, "Producer: gzread error: %s\n", gz_err_msg);
            atomic_store(&global_error, 1);
            free(buffer);
            break;
        }

        Chunk *new_chunk = (Chunk *)malloc(sizeof(Chunk));
        if (new_chunk == NULL) {
            fprintf(stderr, "Producer: Memory allocation failed for chunk\n");
            atomic_store(&global_error, 1);
            free(buffer);
            break;
        }
        new_chunk->buffer = buffer;
        new_chunk->bytes_read = bytes_read;

        pthread_mutex_lock(&queue->mutex);
        ChunkNode *new_node = (ChunkNode *)malloc(sizeof(ChunkNode));
        if (new_node == NULL) {
            fprintf(stderr, "Producer: Memory allocation failed for chunk node\n");
            atomic_store(&global_error, 1);
            free(buffer);
            free(new_chunk);
            pthread_mutex_unlock(&queue->mutex);
            break;
        }
        new_node->chunk = *new_chunk;
        new_node->next = NULL;

        if (queue->tail == NULL) {
            queue->head = new_node;
            queue->tail = new_node;
        } else {
            queue->tail->next = new_node;
            queue->tail = new_node;
        }
        pthread_cond_signal(&queue->cond);
        pthread_mutex_unlock(&queue->mutex);
        
        free(new_chunk); // Free the chunk structure after copying

        if (bytes_read == 0) {
            break; // End of file
        }
    }

    pthread_mutex_lock(&queue->mutex);
    queue->producer_done = 1;
    pthread_cond_broadcast(&queue->cond); // Wake up all consumers
    pthread_mutex_unlock(&queue->mutex);

    return NULL;
}


void *count_newlines(void *arg) {

    ThreadData *data = (ThreadData *)arg;
    ChunkQueue *queue = data->queue;
    size_t local_count = 0;

    while (1) {
        pthread_mutex_lock(&queue->mutex);
        while (queue->head == NULL && !queue->producer_done) {
            pthread_cond_wait(&queue->cond, &queue->mutex);
        }

        if (queue->head == NULL && queue->producer_done) {
            pthread_mutex_unlock(&queue->mutex);
            break; // No more chunks and producer is done
        }

        ChunkNode *node = queue->head;
        queue->head = queue->head->next;
        if (queue->head == NULL) {
            queue->tail = NULL;
        }
        pthread_mutex_unlock(&queue->mutex);

        char *buffer = node->chunk.buffer;
        int bytes_read = node->chunk.bytes_read;

        if (bytes_read == 0) { // End of file chunk
            free(buffer);
            free(node);
            break;
        }

        char *ptr = buffer;
        char *end = buffer + bytes_read;
        while (ptr < end) {
            ptr = (char *)memchr(ptr, '\n', end - ptr);
            if (ptr == NULL) {
                break;
            }
            local_count++;
            ptr++;
        }

        free(buffer);
        free(node);
    }

    data->newline_count = local_count;
    return NULL;
}

int main(int argc, char **argv) {
    if (argc < 2 || argc > 3) {
        fprintf(stderr, "Usage: %s <fastq.gz file> [num_threads]\n", argv[0]);
        return 1;
    }

    int num_threads = NUM_THREADS;
    if (argc == 3) {
        num_threads = atoi(argv[2]);
        if (num_threads <= 0 || num_threads > MAX_THREADS) {
            fprintf(stderr, "Error: Number of threads must be between 1 and %d.\n", MAX_THREADS);
            return 1;
        }
    }


    gzFile file = gzopen(argv[1], "rb");
    if (!file) {
        fprintf(stderr, "Error: Cannot open file %s: %s\n", argv[1], strerror(errno));
        return 1;
    }

    ChunkQueue queue = { .head = NULL, .tail = NULL, .producer_done = 0 };
    pthread_mutex_init(&queue.mutex, NULL);
    pthread_cond_init(&queue.cond, NULL);

    pthread_t producer_thread;
    void *producer_args[2] = {file, &queue};
    if (pthread_create(&producer_thread, NULL, read_chunks, producer_args) != 0) {
        fprintf(stderr, "Error creating producer thread: %s\n", strerror(errno));
        gzclose(file);
        return 1;
    }

    pthread_t consumer_threads[num_threads];
    ThreadData thread_data[num_threads];

    for (int i = 0; i < num_threads; i++) {
        thread_data[i].queue = &queue;
        thread_data[i].newline_count = 0;
        thread_data[i].thread_id = i;

        if (pthread_create(&consumer_threads[i], NULL, count_newlines, &thread_data[i]) != 0) {
            fprintf(stderr, "Error creating consumer thread %d: %s\n", i, strerror(errno));
            
            // Signal producer to stop and join created threads
            atomic_store(&global_error, 1);
            pthread_mutex_lock(&queue.mutex);
            queue.producer_done = 1;
            pthread_cond_broadcast(&queue.cond);
            pthread_mutex_unlock(&queue.mutex);
            
            pthread_join(producer_thread, NULL);
            
            // Join already created consumer threads
            for (int j = 0; j < i; j++) {
                pthread_join(consumer_threads[j], NULL);
            }
            
            gzclose(file);
            pthread_mutex_destroy(&queue.mutex);
            pthread_cond_destroy(&queue.cond);
            return 1;
        }
    }

    pthread_join(producer_thread, NULL);
    for (int i = 0; i < num_threads; i++) {
        pthread_join(consumer_threads[i], NULL);
    }

    gzclose(file);

    pthread_mutex_destroy(&queue.mutex);
    pthread_cond_destroy(&queue.cond);

    if (atomic_load(&global_error)) {
        fprintf(stderr, "An error occurred during processing. Results may be incomplete.\n");
        return 1;
    }

    size_t total_newlines = 0;
    for (int i = 0; i < num_threads; i++) {
        total_newlines += thread_data[i].newline_count;
    }

    size_t sequence_count = total_newlines / 4;
    printf("Total sequences: %zu\n", sequence_count);

    return 0;
}
