#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <zlib.h>
#include <libgen.h>
#include <limits.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/stat.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#define MAX_THREADS 4
#define VERSION "2.0.0"

typedef enum {
    TSV,
    CSV,
    JSON
} output_format_t;

typedef struct {
    char *filepath;
    int abs_path;
    int basename;
    output_format_t output_format;
} task_t;

typedef struct {
    char filepath[PATH_MAX];
    unsigned long total_seqs;
    unsigned long total_len;
    unsigned long n50, n75, n90;
    unsigned long i50;
    double gc_content;
    double avg_len;
    unsigned long min_len, max_len;
} result_t;

pthread_mutex_t io_mutex = PTHREAD_MUTEX_INITIALIZER;
int num_threads = 0;
pthread_mutex_t thread_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t thread_cond = PTHREAD_COND_INITIALIZER;

int compare_desc(const void *a, const void *b) {
    return (*(int *)b - *(int *)a);
}

void *process_file(void *arg) {
    task_t *task = (task_t *)arg;

    gzFile fp;
    if (strcmp(task->filepath, "-") == 0)
        fp = gzdopen(STDIN_FILENO, "r");
    else
        fp = gzopen(task->filepath, "r");
    if (!fp) {
        fprintf(stderr, "Error opening file %s\n", task->filepath);
        pthread_exit(NULL);
    }

    kseq_t *seq = kseq_init(fp);
    unsigned long total_len = 0, total_seqs = 0;
    unsigned long gc_count = 0;
    unsigned long min_len = ULONG_MAX, max_len = 0;
    size_t alloc = 1024;
    unsigned *lengths = malloc(sizeof(unsigned) * alloc);
    if (!lengths) {
        perror("malloc");
        pthread_exit(NULL);
    }

    while (kseq_read(seq) >= 0) {
        unsigned len = seq->seq.l;
        if (total_seqs >= alloc) {
            alloc *= 2;
            unsigned *new_lengths = realloc(lengths, sizeof(unsigned) * alloc);
            if (!new_lengths) {
                perror("realloc");
                free(lengths);
                pthread_exit(NULL);
            }
            lengths = new_lengths;
        }
        lengths[total_seqs++] = len;
        total_len += len;
        if (len < min_len) min_len = len;
        if (len > max_len) max_len = len;
        for (unsigned i = 0; i < len; i++) {
            char c = seq->seq.s[i];
            if (c == 'G' || c == 'g' || c == 'C' || c == 'c') gc_count++;
        }
    }

    kseq_destroy(seq);
    gzclose(fp);

    qsort(lengths, total_seqs, sizeof(unsigned), compare_desc);

    unsigned long sum = 0;
    unsigned long n50 = 0, n75 = 0, n90 = 0, i50 = 0;
    for (unsigned long i = 0; i < total_seqs; i++) {
        sum += lengths[i];
        if (!n50 && sum >= total_len * 0.5) { n50 = lengths[i]; i50 = i + 1; }
        if (!n75 && sum >= total_len * 0.75) n75 = lengths[i];
        if (!n90 && sum >= total_len * 0.90) n90 = lengths[i];
    }

    result_t *res = malloc(sizeof(result_t));
    if (!res) {
        perror("malloc");
        free(lengths);
        pthread_exit(NULL);
    }
    realpath(task->filepath, res->filepath);
    if (task->basename) strcpy(res->filepath, basename(res->filepath));
    res->total_seqs = total_seqs;
    res->total_len = total_len;
    res->n50 = n50;
    res->n75 = n75;
    res->n90 = n90;
    res->i50 = i50;
    res->gc_content = (double)gc_count / total_len * 100.0;
    res->avg_len = (double)total_len / total_seqs;
    res->min_len = min_len;
    res->max_len = max_len;

    free(lengths);

    pthread_mutex_lock(&thread_mutex);
    num_threads--;
    pthread_mutex_unlock(&thread_mutex);

    pthread_exit(res);
}

void print_result(result_t *r, output_format_t fmt) {
    char sep = fmt == CSV ? ',' : '\t';
    printf("%s%c%lu%c%lu%c%lu%c%lu%c%lu%c%lu%c%.2f%c%.2f%c%lu%c%lu\n",
           r->filepath, sep, r->total_seqs, sep, r->total_len, sep, r->n50, sep,
           r->n75, sep, r->n90, sep, r->i50, sep, r->gc_content, sep,
           r->avg_len, sep, r->min_len, sep, r->max_len);
}

void print_json_result(result_t *r, int is_first) {
    if (!is_first) printf(",\n");
    printf("  {\"File\":\"%s\",\"TotSeqs\":%lu,\"TotLen\":%lu,\"N50\":%lu,\"N75\":%lu,\"N90\":%lu,\"I50\":%lu,\"GC\":%.2f,\"Avg\":%.2f,\"Min\":%lu,\"Max\":%lu}",
           r->filepath, r->total_seqs, r->total_len, r->n50, r->n75, r->n90, r->i50, r->gc_content, r->avg_len, r->min_len, r->max_len);
}

void print_help(const char *progname) {
    printf("Usage: %s [options] FILES...\n", progname);
    printf("\nCalculate sequence statistics (N50, GC%%, length stats) for FASTA/FASTQ files.\n\n");
    printf("Arguments:\n");
    printf("  FILES           One or more FASTA/FASTQ files, gzipped or not. '-' for STDIN.\n\n");
    printf("Options:\n");
    printf("  -a, --abs       Print file paths as absolute paths\n");
    printf("  -b, --basename  Print file paths as basename only (e.g., file.fq.gz)\n");
    printf("  -j, --json      Output results in JSON format\n");
    printf("  -c, --csv       Output results in CSV format (default is TSV)\n");
    printf("  -h, --help      Show this help message and exit\n");
    printf("  -v, --version   Show version number and exit\n\n");
    printf("Output Columns (TSV/CSV):\n");
    printf("  Filepath, TotSeqs, TotLen, N50, N75, N90, I50, GC, Avg, Min, Max\n");
    printf("\n");
}

int main(int argc, char *argv[]) {
    int opt, option_index = 0;
    output_format_t output_format = TSV;
    int abs_path = 0, basename_flag = 0;

    static struct option long_opts[] = {
        {"abs", no_argument, 0, 'a'},
        {"basename", no_argument, 0, 'b'},
        {"json", no_argument, 0, 'j'},
        {"csv", no_argument, 0, 'c'},
        {"help", no_argument, 0, 'h'},
        {"version", no_argument, 0, 'v'},
        {0, 0, 0, 0}
    };

    while ((opt = getopt_long(argc, argv, "abjchv", long_opts, &option_index)) != -1) {
        switch (opt) {
            case 'a': abs_path = 1; break;
            case 'b': basename_flag = 1; break;
            case 'j': output_format = JSON; break;
            case 'c': output_format = CSV; break;
            case 'h': print_help(argv[0]); exit(0);
            case 'v': printf("Version: %s\n", VERSION); exit(0);
            default: exit(EXIT_FAILURE);
        }
    }

    int files = argc - optind;
    if (files < 1) {
        fprintf(stderr, "Usage: %s [options] FILES...\n", argv[0]);
        return 1;
    }

    if (output_format == TSV)
        printf("Filepath\tTotSeqs\tTotLen\tN50\tN75\tN90\tI50\tGC\tAvg\tMin\tMax\n");


    // For JSON output, collect all results first
    result_t **all_results = NULL;
    int total_results = 0;
    
    if (output_format == JSON) {
        all_results = malloc(files * sizeof(result_t*));
        if (!all_results) {
            perror("malloc");
            return 1;
        }
    }

    pthread_t threads[MAX_THREADS];
    int running_threads = 0;
    for (int i = optind; i < argc; i++) {
        task_t *t = malloc(sizeof(task_t));
        t->filepath = argv[i];
        t->output_format = output_format;
        t->abs_path = abs_path;
        t->basename = basename_flag;

        pthread_mutex_lock(&thread_mutex);
        while (num_threads >= MAX_THREADS) {
            pthread_mutex_unlock(&thread_mutex);
            usleep(10000);
            pthread_mutex_lock(&thread_mutex);
        }
        pthread_create(&threads[running_threads++], NULL, process_file, t);
        num_threads++;
        pthread_mutex_unlock(&thread_mutex);

        if (running_threads >= MAX_THREADS || i == argc - 1) {
            for (int j = 0; j < running_threads; j++) {
                void *res;
                pthread_join(threads[j], &res);
                if (res) {
                    if (output_format == JSON) {
                        all_results[total_results++] = (result_t *)res;
                    } else {
                        print_result((result_t *)res, output_format);
                        free(res);
                    }
                }
            }
            running_threads = 0;
        }
    }

    // Print JSON array if JSON format was requested
    if (output_format == JSON) {
        printf("[\n");
        for (int i = 0; i < total_results; i++) {
            print_json_result(all_results[i], i == 0);
            free(all_results[i]);
        }
        printf("\n]\n");
        free(all_results);
    }

    return 0;
}

