#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <zlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>
#ifndef CHUNK_SIZE
#define CHUNK_SIZE 1048576 // Default: 1MB chunks
#endif
#define NUM_THREADS 4
typedef struct Chunk {
    char *data;
    size_t size;
    struct Chunk *next;
} Chunk;

typedef struct {
    Chunk *head;
    Chunk *tail;
    int done;
    pthread_mutex_t mutex;
    pthread_cond_t cond;
} ChunkQueue;

typedef struct {
    ChunkQueue *queue;
    size_t seq_count;
    int thread_id;
    size_t chunk_size;
    gzFile file; // only used by producer
} ThreadData;

int global_error = 0;
static void queue_init(ChunkQueue *q) {
    q->head = q->tail = NULL;
    q->done = 0;
    pthread_mutex_init(&q->mutex, NULL);
    pthread_cond_init(&q->cond, NULL);
}

static void queue_destroy(ChunkQueue *q) {
    Chunk *c = q->head;
    while (c) {
        Chunk *next = c->next;
        free(c->data);
        free(c);
        c = next;
    }
    pthread_mutex_destroy(&q->mutex);
    pthread_cond_destroy(&q->cond);
}

static void queue_push(ChunkQueue *q, char *data, size_t size) {
    Chunk *c = malloc(sizeof(Chunk));
    if (!c) {
        fprintf(stderr, "Memory allocation failed in producer\n");
        free(data);
        global_error = 1;
        return;
    }
    c->data = data;
    c->size = size;
    c->next = NULL;

    pthread_mutex_lock(&q->mutex);
    if (q->tail) q->tail->next = c; else q->head = c;
    q->tail = c;
    pthread_cond_signal(&q->cond);
    pthread_mutex_unlock(&q->mutex);
}

static Chunk *queue_pop(ChunkQueue *q) {
    pthread_mutex_lock(&q->mutex);
    while (!q->head && !q->done && !global_error) {
        pthread_cond_wait(&q->cond, &q->mutex);
    }
    Chunk *c = q->head;
    if (c) {
        q->head = c->next;
        if (!q->head) q->tail = NULL;
    }
    pthread_mutex_unlock(&q->mutex);
    return c;
}

static void queue_finish(ChunkQueue *q) {
    pthread_mutex_lock(&q->mutex);
    q->done = 1;
    pthread_cond_broadcast(&q->cond);
    pthread_mutex_unlock(&q->mutex);
}

void *producer(void *arg) {
    ThreadData *data = (ThreadData *)arg;
    while (!global_error) {
        char *buffer = (char *)aligned_alloc(8, data->chunk_size);
        if (!buffer) {
            fprintf(stderr, "Producer: Memory allocation failed\n");
            global_error = 1;
            break;
        }
        int bytes_read = gzread(data->file, buffer, data->chunk_size);
        int gz_err;
        const char *gz_err_msg = gzerror(data->file, &gz_err);
        if (bytes_read < 0 || gz_err != Z_OK) {
            fprintf(stderr, "Producer: gzread error: %s\n", gz_err_msg);
            free(buffer);
            global_error = 1;
            break;
        }
        if (bytes_read == 0) {
            free(buffer);
            break;
        }
        queue_push(data->queue, buffer, (size_t)bytes_read);
    }
    queue_finish(data->queue);
    return NULL;
}

void *consumer(void *arg) {
    ThreadData *data = (ThreadData *)arg;
    size_t local_count = 0;
    while (1) {
        Chunk *c = queue_pop(data->queue);
        if (!c) {
            if (data->queue->done || global_error)
                break;
            else
                continue;
        }
        for (size_t i = 0; i < c->size; i++) {
            if (c->data[i] == '>')
                local_count++;
        }
        free(c->data);
        free(c);
    }
    data->seq_count = local_count;
    return NULL;
}
int main(int argc, char **argv) {
    size_t chunk_size = CHUNK_SIZE;
    static struct option long_opts[] = {
        {"chunk-size", required_argument, 0, 'c'},
        {0, 0, 0, 0}
    };
    int opt;
    while ((opt = getopt_long(argc, argv, "c:", long_opts, NULL)) != -1) {
        switch (opt) {
            case 'c':
                chunk_size = strtoull(optarg, NULL, 10);
                break;
            default:
                fprintf(stderr, "Usage: %s [--chunk-size N] <fasta.gz file>\n", argv[0]);
                return 1;
        }
    }

    if (optind >= argc) {
        fprintf(stderr, "Usage: %s [--chunk-size N] <fasta.gz file>\n", argv[0]);
        return 1;
    }

    gzFile file = gzopen(argv[optind], "rb");
    if (!file) {
        fprintf(stderr, "Error: Cannot open file %s: %s\n", argv[optind], strerror(errno));
        return 1;
    }
    ChunkQueue queue;
    queue_init(&queue);

    pthread_t prod_thread;
    ThreadData prod_data = {.queue = &queue, .chunk_size = chunk_size, .file = file};
    if (pthread_create(&prod_thread, NULL, producer, &prod_data) != 0) {
        fprintf(stderr, "Error creating producer thread: %s\n", strerror(errno));
        gzclose(file);
        queue_destroy(&queue);
        return 1;
    }

    pthread_t threads[NUM_THREADS];
    ThreadData thread_data[NUM_THREADS];
    for (int i = 0; i < NUM_THREADS; i++) {
        thread_data[i].queue = &queue;
        thread_data[i].seq_count = 0;
        thread_data[i].thread_id = i;
        if (pthread_create(&threads[i], NULL, consumer, &thread_data[i]) != 0) {
            fprintf(stderr, "Error creating thread %d: %s\n", i, strerror(errno));
            global_error = 1;
            break;
        }
    }

    pthread_join(prod_thread, NULL);
    for (int i = 0; i < NUM_THREADS; i++) {
        pthread_join(threads[i], NULL);
    }

    queue_destroy(&queue);
    gzclose(file);
    if (global_error) {
        fprintf(stderr, "An error occurred during processing. Results may be incomplete.\n");
        return 1;
    }
    size_t total_sequences = 0;
    for (int i = 0; i < NUM_THREADS; i++) {
        total_sequences += thread_data[i].seq_count;
    }
    printf("Total sequences: %zu\n", total_sequences);
    return 0;
}