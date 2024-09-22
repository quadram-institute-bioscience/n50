#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <pthread.h>
#include <ctype.h>
#include <getopt.h>
#include <stdbool.h>
#include <stdint.h>
#include <limits.h>

// compile time constant DEBUG=0
#define DEBUG 0

#define BUFFER_SIZE 4 * 1024 * 1024  // 4MB buffer
#define MAX_THREADS 16
#define INITIAL_CAPACITY 1000000
#define VERSION "2.0.0"

typedef struct {
    uint64_t *lengths;
    int start;
    int end;
    uint64_t total_length;
} ThreadData;

typedef struct {
    char *filename;
    uint64_t total_length;
    int length_count;
    uint64_t n50;
    uint64_t max_length;
    uint64_t min_length;
    bool is_fastq;
    uint64_t bases[5];  // A, C, G, T, Other
} FileStats;

// Function prototypes
void print_usage(const char *program_name);
void print_version(void);
int compare(const void *a, const void *b);
void *process_chunk(void *arg);
FileStats process_file(const char *filename, bool force_fasta, bool force_fastq, bool extra);

// Main function
int main(int argc, char *argv[]) {
    bool opt_header = false;
    bool opt_n50 = false;
    bool force_fasta = false;
    bool force_fastq = false;
    bool extra = false;
    static struct option long_options[] = {
        {"fasta", no_argument, 0, 'a'},
        {"fastq", no_argument, 0, 'q'},
        {"header", no_argument, 0, 'H'},
        {"n50", no_argument, 0, 'n'},
        {"help", no_argument, 0, 'h'},
        {"version", no_argument, 0, 'v'},
        {"extra", no_argument, 0, 'x'},
        {0, 0, 0, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "aqHnhvx", long_options, NULL)) != -1) {
        switch (opt) {
            case 'a': force_fasta = true; break;
            case 'q': force_fastq = true; break;
            case 'H': opt_header = true; break;
            case 'n': opt_n50 = true; break;
            case 'h': print_usage(argv[0]); return 0;
            case 'v': print_version(); return 0;
            case 'x': extra = true; break;
            default: fprintf(stderr, "Try '%s --help' for more information.\n", argv[0]); return 1;
        }
    }

    if (optind == argc) {
        fprintf(stderr, "Error: No input files specified\n");
        fprintf(stderr, "Try '%s --help' for more information.\n", argv[0]);
        return 1;
    }

    if (opt_header && !opt_n50) {
        if (!extra) {
            printf("Filename\tFormat\tTotal_Length\tTotal_Sequences\tN50\n");
        } else {
            printf("Filename\tFormat\tTotal_Length\tTotal_Sequences\tN50\tMin\tMax\tA\tC\tG\tT\tOther\n");
        }
    }

    for (int i = optind; i < argc; i++) {
        FileStats stats = process_file(argv[i], force_fasta, force_fastq, extra);
        
        if (opt_n50) {
            printf("%s\t%llu\n", stats.filename, stats.n50);
        } else {
            if (extra) {
                double total_bases = stats.bases[0] + stats.bases[1] + stats.bases[2] + stats.bases[3] + stats.bases[4];
                double base_ratios[5];
                for (int j = 0; j < 5; j++) {
                    base_ratios[j] = (double)stats.bases[j] / total_bases * 100;
                }
                printf("%s\t%s\t%llu\t%d\t%llu\t%llu\t%llu\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
                    stats.filename,
                    stats.is_fastq ? "FASTQ" : "FASTA",
                    stats.total_length,
                    stats.length_count,
                    stats.n50,
                    stats.min_length,
                    stats.max_length,
                    base_ratios[0],
                    base_ratios[1],
                    base_ratios[2],
                    base_ratios[3],
                    base_ratios[4]);
            } else {
                printf("%s\t%s\t%llu\t%d\t%llu\n", stats.filename, stats.is_fastq ? "FASTQ" : "FASTA", stats.total_length, stats.length_count, stats.n50);
            }
        }

        free(stats.filename);
    }

    return 0;
}

void print_usage(const char *program_name) {
    printf("Usage: %s [OPTIONS] [FILENAME...]\n", program_name);
    printf("\nOptions:\n");
    printf("  -a, --fasta        Force FASTA input format\n");
    printf("  -q, --fastq        Force FASTQ input format\n");
    printf("  -H, --header       Print header in output\n");
    printf("  -n, --n50          Output only N50 value\n");
    printf("  -x, --extra        Output extra statistics\n");
    printf("  -h, --help         Display this help message and exit\n");
    printf("      --version      Display version information and exit\n");
    printf("\nDescription:\n");
    printf("  Calculate N50 and other sequence statistics from FASTA or FASTQ files.\n");
    printf("  Supports multiple input files, automatic format detection, and reading from STDIN.\n");
    printf("  Use '-' as filename to read uncompressed input from STDIN.\n");
}

void print_version(void) {
    printf("N50 Calculator version %s\n", VERSION);
    printf("Copyright (C) 2024 Andrea Telatin\n");
    printf("License: MIT\n");
}

int compare(const void *a, const void *b) {
    return (*(uint64_t*)b - *(uint64_t*)a);
}

void *process_chunk(void *arg) {
    ThreadData *data = (ThreadData*)arg;
    uint64_t local_total_length = 0;
    uint64_t *local_lengths = malloc((data->end - data->start) * sizeof(uint64_t));
    int local_count = 0;

    for (int i = data->start; i < data->end; i++) {
        local_total_length += data->lengths[i];
        local_lengths[local_count++] = data->lengths[i];
    }

    ThreadData *result = malloc(sizeof(ThreadData));
    result->lengths = local_lengths;
    result->start = data->start;
    result->end = data->end;
    result->total_length = local_total_length;

    return result;
}


FileStats process_file(const char *filename, bool force_fasta, bool force_fastq, bool extra) {
    FILE *fp = NULL;
    gzFile gzfp = NULL;
    bool is_stdin = (strcmp(filename, "-") == 0);
    bool is_gzipped = (!is_stdin && strstr(filename, ".gz") != NULL);

    if (is_stdin) {
        fp = stdin;
    } else if (is_gzipped) {
        gzfp = gzopen(filename, "r");
        if (!gzfp) {
            fprintf(stderr, "Error: Cannot open file %s\n", filename);
            exit(1);
        }
    } else {
        fp = fopen(filename, "r");
        if (!fp) {
            fprintf(stderr, "Error: Cannot open file %s\n", filename);
            exit(1);
        }
    }

    char buffer[BUFFER_SIZE];
    uint64_t *chunk_lengths = NULL;
    int chunk_count = 0;
    int chunk_capacity = INITIAL_CAPACITY;

    // Extra statistics
    uint64_t current_length = 0;
    uint64_t min = UINT64_MAX;
    uint64_t max = 0;
    uint64_t bases[5] = {0};  // A, C, G, T, Other

    bool is_fastq = force_fastq || (!force_fasta && (strstr(filename, ".fastq") != NULL || strstr(filename, ".fq") != NULL));

    chunk_lengths = malloc(chunk_capacity * sizeof(uint64_t));

    // lookup table
        // Create a lookup table for base counting
    uint8_t base_lookup[256] = {0};
    base_lookup['A'] = base_lookup['a'] = 0;
    base_lookup['C'] = base_lookup['c'] = 1;
    base_lookup['G'] = base_lookup['g'] = 2;
    base_lookup['T'] = base_lookup['t'] = 3;

    if (is_fastq) {

        #if DEBUG
            fprintf(stderr, "Debug: Processing FASTQ file %s\n", filename);
        #endif
        int newline_count = 0;
        bool in_sequence = false;
        while (1) {
            size_t bytes_read;
            if (is_gzipped) {
                bytes_read = gzread(gzfp, buffer, BUFFER_SIZE);
            } else {
                bytes_read = fread(buffer, 1, BUFFER_SIZE, fp);
            }
            if (bytes_read == 0) break;

            for (size_t i = 0; i < bytes_read; i++) {
                char c = buffer[i];
                if (c == '\n') {
                    newline_count++;
                    if (newline_count % 4 == 1) {
                        in_sequence = true;
                    } else if (newline_count % 4 == 2) {
                        in_sequence = false;
                        if (chunk_count == chunk_capacity) {
                            chunk_capacity *= 2;
                            chunk_lengths = realloc(chunk_lengths, chunk_capacity * sizeof(uint64_t));
                        }
                        chunk_lengths[chunk_count++] = current_length;
                        
                        if (extra) {
                            if (current_length > max) max = current_length;
                            if (current_length < min) min = current_length;
                        }
                        current_length = 0;
                    }
                } else if (in_sequence) {
                    current_length++;
                    if (extra) {
                        bases[base_lookup[(uint8_t)c]]++;
                    }
                }
            }
        }
    } else {  // FASTA
        bool in_sequence = false;
        while (1) {
            char *result;
            if (is_gzipped) {
                result = gzgets(gzfp, buffer, BUFFER_SIZE);
            } else {
                result = fgets(buffer, BUFFER_SIZE, fp);
            }
            if (!result) break;

            if (buffer[0] == '>') {
                if (current_length > 0) {
                    if (chunk_count == chunk_capacity) {
                        chunk_capacity *= 2;
                        chunk_lengths = realloc(chunk_lengths, chunk_capacity * sizeof(uint64_t));
                    }
                    chunk_lengths[chunk_count++] = current_length;
                    if (extra) {
                        if (current_length > max) max = current_length;
                        if (current_length < min) min = current_length;
                    }
                    current_length = 0;
                }
                in_sequence = true;
            } else if (in_sequence) {
                uint64_t len = strcspn(buffer, "\n");
                current_length += len;
                if (extra) {
                    for (uint64_t i = 0; i < len; i++) {
                        bases[base_lookup[(uint8_t)buffer[i]]]++;
                    }
                }
            }
        }
        if (current_length > 0) {
            if (chunk_count == chunk_capacity) {
                chunk_capacity *= 2;
                chunk_lengths = realloc(chunk_lengths, chunk_capacity * sizeof(uint64_t));
            }
            chunk_lengths[chunk_count++] = current_length;
            if (extra) {
                if (current_length > max) max = current_length;
                if (current_length < min) min = current_length;
            }
        }
    }

    if (is_gzipped) {
        gzclose(gzfp);
    } else if (!is_stdin) {
        fclose(fp);
    }

   // Calculate total length and other statistics
    uint64_t total_length = 0;
    for (int i = 0; i < chunk_count; i++) {
        total_length += chunk_lengths[i];
    }

    // Calculate N50
    qsort(chunk_lengths, chunk_count, sizeof(uint64_t), compare);
    uint64_t cumulative_length = 0;
    uint64_t n50 = 0;
    for (int i = chunk_count - 1; i >= 0; i--) {
        cumulative_length += chunk_lengths[i];
        if (cumulative_length >= total_length / 2) {
            n50 = chunk_lengths[i];
            break;
        }
    }

    free(chunk_lengths);

    FileStats stats = {
        .filename = strdup(filename),
        .total_length = total_length,
        .length_count = chunk_count,
        .n50 = n50,
        .is_fastq = is_fastq,
        .max_length = max,
        .min_length = min
    };
    
    if (extra) {
        for (int i = 0; i < 5; i++) {
            stats.bases[i] = bases[i];
        }
    }

    #if DEBUG
        fprintf(stderr, "Debug: Final stats - total_length: %lu, chunk_count: %d\n", stats.total_length, stats.length_count);
    #endif
    return stats;
}