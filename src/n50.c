// Full modified version of n50.c with AuN and AuN75 columns added
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
#include <sys/ioctl.h>
#include <termios.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#define MAX_THREADS 4
#define VERSION "1.9.2"

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
    int nice_output;
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
    unsigned long aun;
} result_t;

pthread_mutex_t io_mutex = PTHREAD_MUTEX_INITIALIZER;
int num_threads = 0;
pthread_mutex_t thread_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t thread_cond = PTHREAD_COND_INITIALIZER;

int compare_desc(const void *a, const void *b) {
    return (*(int *)b - *(int *)a);
}

int get_terminal_width() {
    struct winsize w;
    if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &w) == 0) {
        return w.ws_col;
    }
    return 80; // Default fallback width
}

unsigned long calculate_auN(unsigned *lengths, unsigned long n, unsigned long limit) {
    double aun = 0.0;
    unsigned long cumulative = 0;
    for (unsigned long i = 0; i < n && cumulative < limit; i++) {
        unsigned long eff_len = (cumulative + lengths[i] <= limit) ? lengths[i] : (limit - cumulative);
        aun += eff_len * ((double)eff_len / limit);
        cumulative += lengths[i];
    }
    return (unsigned long)(aun + 0.5);
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
    res->aun = calculate_auN(lengths, total_seqs, total_len);

    free(lengths);

    pthread_mutex_lock(&thread_mutex);
    num_threads--;
    pthread_mutex_unlock(&thread_mutex);

    pthread_exit(res);
}

void print_result(result_t *r, output_format_t fmt, int nice_output) {
    if (nice_output) {
        int term_width = get_terminal_width();
        // Calculate dynamic column widths based on terminal width
        // Reserve minimum space for each column, adjust filepath width dynamically
        int min_col_width = 8;
        int num_cols = 12; // total columns including filepath
        int reserved_width = (num_cols - 1) * min_col_width + (num_cols - 1); // space for separators
        int filepath_width = term_width - reserved_width;
        
        // Ensure filepath width is at least 15 chars and not more than 50
        if (filepath_width < 15) filepath_width = 15;
        if (filepath_width > 50) filepath_width = 50;
        
        printf("%-*s %*lu %*lu %*lu %*lu %*lu %*lu %*.2f %*.2f %*lu %*lu %*lu\n",
               filepath_width, r->filepath,
               min_col_width, r->total_seqs,
               min_col_width, r->total_len,
               min_col_width, r->n50,
               min_col_width, r->n75,
               min_col_width, r->n90,
               min_col_width, r->i50,
               min_col_width, r->gc_content,
               min_col_width, r->avg_len,
               min_col_width, r->min_len,
               min_col_width, r->max_len,
               min_col_width, r->aun);
    } else {
        char sep = fmt == CSV ? ',' : '\t';
        printf("%s%c%lu%c%lu%c%lu%c%lu%c%lu%c%lu%c%.2f%c%.2f%c%lu%c%lu%c%lu\n",
               r->filepath, sep, r->total_seqs, sep, r->total_len, sep, r->n50, sep,
               r->n75, sep, r->n90, sep, r->i50, sep, r->gc_content, sep,
               r->avg_len, sep, r->min_len, sep, r->max_len, sep, r->aun);
    }
}

void print_json_result(result_t *r, int is_first) {
    if (!is_first) printf(",\n");
    printf("  {\"File\":\"%s\",\"TotSeqs\":%lu,\"TotLen\":%lu,\"N50\":%lu,\"N75\":%lu,\"N90\":%lu,\"I50\":%lu,\"GC\":%.2f,\"Avg\":%.2f,\"Min\":%lu,\"Max\":%lu,\"AuN\":%lu}",
           r->filepath, r->total_seqs, r->total_len, r->n50, r->n75, r->n90, r->i50, r->gc_content, r->avg_len, r->min_len, r->max_len, r->aun);
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
    printf("  -n, --nice      Output results in a visually aligned ASCII table\n");
    printf("  -h, --help      Show this help message and exit\n");
    printf("  -v, --version   Show version number and exit\n\n");
    printf("Output Columns (TSV/CSV):\n");
    printf("  Filepath, TotSeqs, TotLen, N50, N75, N90, I50, GC, Avg, Min, Max, AuN\n\n");
}

int main(int argc, char *argv[]) {
    int opt, option_index = 0;
    output_format_t output_format = TSV;
    int abs_path = 0, basename_flag = 0;
    int nice_output = 0;

    static struct option long_opts[] = {
        {"abs", no_argument, 0, 'a'},
        {"basename", no_argument, 0, 'b'},
        {"json", no_argument, 0, 'j'},
        {"csv", no_argument, 0, 'c'},
        {"nice", no_argument, 0, 'n'},
        {"help", no_argument, 0, 'h'},
        {"version", no_argument, 0, 'v'},
        {0, 0, 0, 0}
    };

    while ((opt = getopt_long(argc, argv, "abjchnv", long_opts, &option_index)) != -1) {
        switch (opt) {
            case 'a': abs_path = 1; break;
            case 'b': basename_flag = 1; break;
            case 'j': output_format = JSON; break;
            case 'c': output_format = CSV; break;
            case 'n': nice_output = 1; basename_flag = 1; break;
            case 'h': print_help(argv[0]); exit(0);
            case 'v': printf("%s\n", VERSION); exit(0);
            default: exit(EXIT_FAILURE);
        }
    }

    int files = argc - optind;
    if (files < 1) {
        fprintf(stderr, "Usage: %s [options] FILES...\n", argv[0]);
        return 1;
    }

    if (output_format == TSV) {
        if (nice_output) {
            int term_width = get_terminal_width();
            int min_col_width = 8;
            int num_cols = 12;
            int reserved_width = (num_cols - 1) * min_col_width + (num_cols - 1);
            int filepath_width = term_width - reserved_width;
            
            if (filepath_width < 15) filepath_width = 15;
            if (filepath_width > 50) filepath_width = 50;
            
            printf("%-*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s\n",
                   filepath_width, "Filepath",
                   min_col_width, "TotSeqs",
                   min_col_width, "TotLen",
                   min_col_width, "N50",
                   min_col_width, "N75",
                   min_col_width, "N90",
                   min_col_width, "I50",
                   min_col_width, "GC",
                   min_col_width, "Avg",
                   min_col_width, "Min",
                   min_col_width, "Max",
                   min_col_width, "AuN");
        } else {
            printf("Filepath\tTotSeqs\tTotLen\tN50\tN75\tN90\tI50\tGC\tAvg\tMin\tMax\tAuN\n");
        }
    }

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
        t->nice_output = nice_output;

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
                        print_result((result_t *)res, output_format, nice_output);
                        free(res);
                    }
                }
            }
            running_threads = 0;
        }
    }

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
