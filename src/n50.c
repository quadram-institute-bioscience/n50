#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <pthread.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

#define BUF_SIZE (100 * 1024 * 1024)
#define MAX_PATH 4096
#define MAX_SEQ_LEN (10 * 1024 * 1024)

typedef struct {
    const char *filename;
    unsigned long tot_seqs;
    unsigned long tot_bp;
    unsigned *lengths;
    unsigned len_cap;
    unsigned len_size;
} FileStat;

void *process_file(void *arg);
int compare_desc(const void *a, const void *b);
int is_fastq(const char *filename);

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s FILES...\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    pthread_t threads[argc - 1];
    FileStat stats[argc - 1];

    for (int i = 1; i < argc; ++i) {
        stats[i - 1].filename = argv[i];
        stats[i - 1].tot_seqs = 0;
        stats[i - 1].tot_bp = 0;
        stats[i - 1].lengths = malloc(sizeof(unsigned) * 1000000);
        stats[i - 1].len_cap = 1000000;
        stats[i - 1].len_size = 0;
        pthread_create(&threads[i - 1], NULL, process_file, &stats[i - 1]);
    }

    printf("Filename\tTotSeqs\tTotBp\tN50\n");
    for (int i = 0; i < argc - 1; ++i) {
        pthread_join(threads[i], NULL);
        qsort(stats[i].lengths, stats[i].len_size, sizeof(unsigned), compare_desc);
        unsigned long sum = 0, half = stats[i].tot_bp / 2;
        unsigned N50 = 0;
        for (unsigned j = 0; j < stats[i].len_size; ++j) {
            sum += stats[i].lengths[j];
            if (sum >= half) {
                N50 = stats[i].lengths[j];
                break;
            }
        }
        printf("%s\t%lu\t%lu\t%u\n", stats[i].filename, stats[i].tot_seqs, stats[i].tot_bp, N50);
        free(stats[i].lengths);
    }
    return 0;
}

void *process_file(void *arg) {
    FileStat *stat = (FileStat *)arg;
    int fastq = is_fastq(stat->filename);
    gzFile fp = gzopen(stat->filename, "rb");
    if (!fp) {
        perror(stat->filename);
        pthread_exit(NULL);
    }
    char *buf = malloc(BUF_SIZE);
    char *seq = malloc(MAX_SEQ_LEN);
    size_t seq_len = 0;
    int state = fastq ? 0 : 1;  // FASTQ state machine: 0-1-2-3
    while (gzgets(fp, buf, BUF_SIZE)) {
        size_t len = strlen(buf);
        if (buf[len - 1] == '\n') buf[len - 1] = '\0';
        if (fastq) {
            if (state == 1) {
                seq_len = strlen(buf);
                if (stat->len_size >= stat->len_cap) {
                    stat->len_cap *= 2;
                    stat->lengths = realloc(stat->lengths, stat->len_cap * sizeof(unsigned));
                }
                stat->lengths[stat->len_size++] = seq_len;
                stat->tot_bp += seq_len;
                stat->tot_seqs++;
            }
            state = (state + 1) % 4;
        } else {
            if (buf[0] == '>') {
                if (seq_len > 0) {
                    if (stat->len_size >= stat->len_cap) {
                        stat->len_cap *= 2;
                        stat->lengths = realloc(stat->lengths, stat->len_cap * sizeof(unsigned));
                    }
                    stat->lengths[stat->len_size++] = seq_len;
                    stat->tot_bp += seq_len;
                    stat->tot_seqs++;
                    seq_len = 0;
                }
            } else {
                seq_len += strlen(buf);
            }
        }
    }
    if (!fastq && seq_len > 0) {
        if (stat->len_size >= stat->len_cap) {
            stat->len_cap *= 2;
            stat->lengths = realloc(stat->lengths, stat->len_cap * sizeof(unsigned));
        }
        stat->lengths[stat->len_size++] = seq_len;
        stat->tot_bp += seq_len;
        stat->tot_seqs++;
    }
    gzclose(fp);
    free(buf);
    free(seq);
    pthread_exit(NULL);
}

int compare_desc(const void *a, const void *b) {
    return (*(unsigned *)b) - (*(unsigned *)a);
}

int is_fastq(const char *filename) {
    const char *dot = strrchr(filename, '.');
    if (!dot) return 0;
    if (strcmp(dot, ".gz") == 0) {
        char *dup = strdup(filename);
        char *dot2 = strrchr(dup, '.');
        if (dot2) *dot2 = '\0';
        dot2 = strrchr(dup, '.');
        int isfq = (dot2 && (!strcmp(dot2, ".fq") || !strcmp(dot2, ".fastq")));
        free(dup);
        return isfq;
    }
    return (!strcmp(dot, ".fq") || !strcmp(dot, ".fastq"));
}

