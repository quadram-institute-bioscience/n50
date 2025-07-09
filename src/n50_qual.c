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
#include <math.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#define MAX_THREADS 4
#define VERSION "1.9.4"

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
    int qual_offset;
    char *output_file;
} task_t;

typedef struct {
    char *readname;
    unsigned int length;
    double avg_quality;
} seq_qual_t;

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
    unsigned long total_quality;
    unsigned long q20_count;
    unsigned long q30_count;
    double avg_quality;
    double q20_fraction;
    double q30_fraction;
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

void write_seq_qual_tsv(seq_qual_t *seq_quals, unsigned long total_seqs, const char *output_file) {
    FILE *fp = fopen(output_file, "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open output file %s\n", output_file);
        return;
    }
    
    // Write header
    fprintf(fp, "readname\tlength\tavg_qual\n");
    
    // Write data
    for (unsigned long i = 0; i < total_seqs; i++) {
        fprintf(fp, "%s\t%u\t%.2f\n", 
                seq_quals[i].readname ? seq_quals[i].readname : "unknown",
                seq_quals[i].length,
                seq_quals[i].avg_quality);
    }
    
    fclose(fp);
}

void free_seq_quals(seq_qual_t *seq_quals, unsigned long total_seqs) {
    for (unsigned long i = 0; i < total_seqs; i++) {
        if (seq_quals[i].readname) {
            free(seq_quals[i].readname);
        }
    }
    free(seq_quals);
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
    int first_seq = 1;
    unsigned long total_len = 0, total_seqs = 0;
    unsigned long gc_count = 0;
    unsigned long min_len = ULONG_MAX, max_len = 0;
    unsigned long total_quality = 0, q20_count = 0, q30_count = 0;
    double total_error_prob_sum = 0.0;
    size_t alloc = 1024;
    unsigned *lengths = malloc(sizeof(unsigned) * alloc);
    seq_qual_t *seq_quals = malloc(sizeof(seq_qual_t) * alloc);
    if (!lengths || !seq_quals) {
        perror("malloc");
        if (lengths) free(lengths);
        if (seq_quals) free_seq_quals(seq_quals, 0);
        pthread_exit(NULL);
    }

    while (kseq_read(seq) >= 0) {
        // Check if this is FASTA format (no quality scores)
        if (first_seq && seq->qual.l == 0) {
            fprintf(stderr, "Error: File %s appears to be in FASTA format. This tool requires FASTQ files with quality scores.\n", task->filepath);
            kseq_destroy(seq);
            gzclose(fp);
            free(lengths);
            free_seq_quals(seq_quals, 0);
            pthread_exit(NULL);
        }
        first_seq = 0;
        
        unsigned len = seq->seq.l;
        if (total_seqs >= alloc) {
            alloc *= 2;
            unsigned *new_lengths = realloc(lengths, sizeof(unsigned) * alloc);
            seq_qual_t *new_seq_quals = realloc(seq_quals, sizeof(seq_qual_t) * alloc);
            if (!new_lengths || !new_seq_quals) {
                perror("realloc");
                free(lengths);
                free_seq_quals(seq_quals, total_seqs);
                pthread_exit(NULL);
            }
            lengths = new_lengths;
            seq_quals = new_seq_quals;
        }
        lengths[total_seqs] = len;
        total_len += len;
        if (len < min_len) min_len = len;
        if (len > max_len) max_len = len;
        
        // Parse quality values
        double seq_error_prob_sum = 0.0;
        for (unsigned i = 0; i < len; i++) {
            char c = seq->seq.s[i];
            if (c == 'G' || c == 'g' || c == 'C' || c == 'c') gc_count++;
            
            // Parse quality score
            int qual = (int)seq->qual.s[i] - task->qual_offset;
            // Convert quality to error probability: P = 10^(-Q/10)
            double error_prob = pow(10.0, -qual / 10.0);
            seq_error_prob_sum += error_prob;
            total_error_prob_sum += error_prob;
            total_quality += qual;
            if (qual >= 20) q20_count++;
            if (qual >= 30) q30_count++;
        }
        
        // Store sequence length, average quality, and readname
        seq_quals[total_seqs].length = len;
        // Calculate average quality using logarithmic method: Q_avg = -10 * log10(P_avg)
        double avg_error_prob = seq_error_prob_sum / len;
        seq_quals[total_seqs].avg_quality = (avg_error_prob == 0.0) ? 0.0 : -10.0 * log10(avg_error_prob);
        seq_quals[total_seqs].readname = malloc(strlen(seq->name.s) + 1);
        if (seq_quals[total_seqs].readname) {
            strcpy(seq_quals[total_seqs].readname, seq->name.s);
        }
        
        total_seqs++;
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
        free_seq_quals(seq_quals, total_seqs);
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
    res->total_quality = total_quality;
    res->q20_count = q20_count;
    res->q30_count = q30_count;
    // Calculate average quality using logarithmic method: Q_avg = -10 * log10(P_avg)
    double avg_error_prob = total_error_prob_sum / total_len;
    res->avg_quality = (avg_error_prob == 0.0) ? 0.0 : -10.0 * log10(avg_error_prob);
    res->q20_fraction = (double)q20_count / total_len;
    res->q30_fraction = (double)q30_count / total_len;

    // Write TSV output if requested
    if (task->output_file) {
        write_seq_qual_tsv(seq_quals, total_seqs, task->output_file);
    }

    free(lengths);
    free_seq_quals(seq_quals, total_seqs);

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
        
        printf("%-*s %*lu %*lu %*lu %*lu %*lu %*lu %*.2f %*.2f %*lu %*lu %*lu %*.2f %*.2f %*.2f\n",
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
               min_col_width, r->aun,
               min_col_width, r->avg_quality,
               min_col_width, r->q20_fraction * 100.0,
               min_col_width, r->q30_fraction * 100.0);
    } else {
        char sep = fmt == CSV ? ',' : '\t';
        printf("%s%c%lu%c%lu%c%lu%c%lu%c%lu%c%lu%c%.2f%c%.2f%c%lu%c%lu%c%lu%c%.2f%c%.2f%c%.2f\n",
               r->filepath, sep, r->total_seqs, sep, r->total_len, sep, r->n50, sep,
               r->n75, sep, r->n90, sep, r->i50, sep, r->gc_content, sep,
               r->avg_len, sep, r->min_len, sep, r->max_len, sep, r->aun, sep,
               r->avg_quality, sep, r->q20_fraction * 100.0, sep, r->q30_fraction * 100.0);
    }
}

void print_json_result(result_t *r, int is_first) {
    if (!is_first) printf(",\n");
    printf("  {\"File\":\"%s\",\"TotSeqs\":%lu,\"TotLen\":%lu,\"N50\":%lu,\"N75\":%lu,\"N90\":%lu,\"I50\":%lu,\"GC\":%.2f,\"Avg\":%.2f,\"Min\":%lu,\"Max\":%lu,\"AuN\":%lu,\"AvgQual\":%.2f,\"Q20\":%.2f,\"Q30\":%.2f}",
           r->filepath, r->total_seqs, r->total_len, r->n50, r->n75, r->n90, r->i50, r->gc_content, r->avg_len, r->min_len, r->max_len, r->aun, r->avg_quality, r->q20_fraction * 100.0, r->q30_fraction * 100.0);
}

void print_help(const char *progname) {
    printf("Usage: %s [options] FILES...\n", progname);
    printf("\nCalculate sequence and quality statistics for FASTQ files.\n\n");
    printf("Arguments:\n");
    printf("  FILES           One or more FASTQ files, gzipped or not. '-' for STDIN.\n\n");
    printf("Options:\n");
    printf("  -a, --abs       Print file paths as absolute paths\n");
    printf("  -b, --basename  Print file paths as basename only (e.g., file.fq.gz)\n");
    printf("  -j, --json      Output results in JSON format\n");
    printf("  -c, --csv       Output results in CSV format (default is TSV)\n");
    printf("  -n, --nice      Output results in a visually aligned ASCII table\n");
    printf("  -o, --output FILE  Save per-sequence data (readname, length, avg_qual) to TSV file\n");
    printf("  --offset INT    Phred quality score offset (default: 33)\n");
    printf("  -h, --help      Show this help message and exit\n");
    printf("  -v, --version   Show version number and exit\n\n");
    printf("Output Columns (TSV/CSV):\n");
    printf("  Filepath, TotSeqs, TotLen, N50, N75, N90, I50, GC, Avg, Min, Max, AuN, AvgQual, Q20, Q30\n\n");
}

int main(int argc, char *argv[]) {
    int opt, option_index = 0;
    output_format_t output_format = TSV;
    int abs_path = 0, basename_flag = 0;
    int nice_output = 0;
    int qual_offset = 33;
    char *output_file = NULL;

    static struct option long_opts[] = {
        {"abs", no_argument, 0, 'a'},
        {"basename", no_argument, 0, 'b'},
        {"json", no_argument, 0, 'j'},
        {"csv", no_argument, 0, 'c'},
        {"nice", no_argument, 0, 'n'},
        {"output", required_argument, 0, 'o'},
        {"offset", required_argument, 0, 'O'},
        {"help", no_argument, 0, 'h'},
        {"version", no_argument, 0, 'v'},
        {0, 0, 0, 0}
    };

    while ((opt = getopt_long(argc, argv, "abjchno:O:v", long_opts, &option_index)) != -1) {
        switch (opt) {
            case 'a': abs_path = 1; break;
            case 'b': basename_flag = 1; break;
            case 'j': output_format = JSON; break;
            case 'c': output_format = CSV; break;
            case 'n': nice_output = 1; basename_flag = 1; break;
            case 'o': output_file = optarg; break;
            case 'O': qual_offset = atoi(optarg); break;
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
            
            printf("%-*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s\n",
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
                   min_col_width, "AuN",
                   min_col_width, "AvgQual",
                   min_col_width, "Q20",
                   min_col_width, "Q30");
        } else {
            printf("Filepath\tTotSeqs\tTotLen\tN50\tN75\tN90\tI50\tGC\tAvg\tMin\tMax\tAuN\tAvgQual\tQ20\tQ30\n");
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
        t->qual_offset = qual_offset;
        t->output_file = output_file;

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
