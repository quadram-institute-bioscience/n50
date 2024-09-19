#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
/*
 * DNA Sequence Generator
 * Andrea Telatin 2023, (C) Quadram Institute Bioscience
 * 
 * This program generates random DNA sequences in FASTA or FASTQ format.
 * It creates multiple files containing sequences of varying lengths,
 * simulating genomic data for testing purposes.
 *
 * Features:
 * - Generates sequences with customizable length ranges and sequence counts
 * - Calculates and includes N50 statistics in the output filenames
 * - Supports both FASTA and FASTQ output formats
 * - Output filename format: N50_TOTSEQS_SUMLEN.{fasta|fastq}
 *
 * Usage: ./program <min_seqs> <max_seqs> <min_len> <max_len> <tot_files> <format> <outdir>
 * 
 * This tool is useful for testing sequence analysis software, benchmarking
 * bioinformatics pipelines, and generating sample data for educational purposes.
 */
#define MAX_SEQ_NAME_LEN 256
#define MAX_FILENAME_LEN 256

// Function prototypes
long long calculate_n50(const int *lengths, int num_seqs, long long *total_length);
void write_fasta(char **sequences, int *lengths, int num_seqs, const char *outfile);
void write_fastq(char **sequences, int *lengths, int num_seqs, const char *outfile);

void gen_ctg_len(int min_len, int max_len, int num_seqs, int *lengths);
int *generate_contigs(int N50, int SUM_LEN, int TOT_SEQS);
void generate_sequence(int length, char *sequence);
void free_resources(char **sequences, int *lengths, int num_seqs);

int main(int argc, char *argv[]) {
    if (argc != 8) {
        fprintf(stderr, "Usage: %s <min_seqs> <max_seqs> <min_len> <max_len> <tot_files> <format> <outdir>\n", argv[0]);
        return 1;
    }

    int min_seqs = atoi(argv[1]);
    int max_seqs = atoi(argv[2]);
    int min_len = atoi(argv[3]);
    int max_len = atoi(argv[4]);
    int tot_files = atoi(argv[5]);
    char *format = argv[6];
    char *outdir = argv[7];

    // Input validation
    if (min_seqs <= 0 || max_seqs <= 0 || min_len <= 0 || max_len <= 0 || tot_files <= 0) {
        fprintf(stderr, "Error: All numeric inputs must be positive integers.\n");
        return 1;
    }

    if (min_seqs > max_seqs || min_len > max_len) {
        fprintf(stderr, "Error: Min values must be less than or equal to max values.\n");
        return 1;
    }

    // Convert format to lowercase for case-insensitive comparison
    for (char *p = format; *p; ++p) *p = tolower(*p);
    if (strcmp(format, "fasta") != 0 && strcmp(format, "fastq") != 0) {
        fprintf(stderr, "Error: Invalid format. Use 'fasta' or 'fastq'.\n");
        return 1;
    }

    // Log input parameters
    fprintf(stderr, "Parameters: min_seqs=%d, max_seqs=%d, min_len=%d, max_len=%d, tot_files=%d, format=%s, outdir=%s\n",
            min_seqs, max_seqs, min_len, max_len, tot_files, format, outdir);

    int seed = 42;  // Fixed seed for reproducibility
    srand(seed);

    for (int i = 1; i <= tot_files; i++) {
        int num_seqs = rand() % (max_seqs - min_seqs + 1) + min_seqs;
        int *contig_lengths = malloc(num_seqs * sizeof(int));
        if (!contig_lengths) {
            fprintf(stderr, "Memory allocation failed for contig_lengths\n");
            return 1;
        }
        fprintf(stderr, "%d file [%d num seqs]\n", i, num_seqs);
        gen_ctg_len(min_len, max_len, num_seqs, contig_lengths);
        
        long long total_length; // Change to long long
        int N50 = calculate_n50(contig_lengths, num_seqs, &total_length);
        fprintf(stderr, "\tTotal length: %lld\n", total_length); // Use %lld for long long
        fprintf(stderr, "\tN50: %d\n", N50);
        char outfile[MAX_FILENAME_LEN];
        snprintf(outfile, sizeof(outfile), "%s/%d_%d_%lld.%s", outdir, N50, num_seqs, total_length, format); // Use %lld for long long

        char **sequences = malloc(num_seqs * sizeof(char *));
        if (!sequences) {
            fprintf(stderr, "Memory allocation failed for sequences\n");
            free(contig_lengths);
            return 1;
        }

        for (int j = 0; j < num_seqs; j++) {
            sequences[j] = malloc((contig_lengths[j] + 1) * sizeof(char));
            if (!sequences[j]) {
                fprintf(stderr, "Memory allocation failed for sequence %d\n", j);
                free_resources(sequences, contig_lengths, j);
                return 1;
            }
            generate_sequence(contig_lengths[j], sequences[j]);
        }

        if (strcmp(format, "fasta") == 0) {
            write_fasta(sequences, contig_lengths, num_seqs, outfile);
        } else {
            write_fastq(sequences, contig_lengths, num_seqs, outfile);
        }
        fprintf(stderr, "  [Done: %s]\n", outfile);
        free_resources(sequences, contig_lengths, num_seqs);
    }

    return 0;
}

void gen_ctg_len(int min_len, int max_len, int num_seqs, int *lengths) {
    for (int i = 0; i < num_seqs; i++) {
        lengths[i] = rand() % (max_len - min_len + 1) + min_len;
    }
}

int compare_ints(const void *a, const void *b) {
    return (*(int*)b - *(int*)a);
}
long long calculate_n50(const int *lengths, int num_seqs, long long *total_length) {
    int sorted_lengths[num_seqs];
    memcpy(sorted_lengths, lengths, sizeof(int) * num_seqs);
    
    qsort(sorted_lengths, num_seqs, sizeof(int), compare_ints);
    
    *total_length = 0;
    for (int i = 0; i < num_seqs; i++) {
        *total_length += sorted_lengths[i];
    }
    
    long long cumulative_length = 0;
    for (int i = 0; i < num_seqs; i++) {
        cumulative_length += sorted_lengths[i];
        if (cumulative_length >= *total_length / 2) {
            return sorted_lengths[i];
        }
    }
    
    return -1;  // Should never reach this point
}

int *generate_contigs(int N50, int SUM_LEN, int TOT_SEQS) {
    if (N50 > SUM_LEN || TOT_SEQS < 1) {
        fprintf(stderr, "Invalid input: N50 must be <= SUM_LEN and TOT_SEQS >= 1\n");
        return NULL;
    }
    
    int *contig_list = malloc(TOT_SEQS * sizeof(int));
    if (!contig_list) {
        fprintf(stderr, "Memory allocation failed\n");
        return NULL;
    }
    
    int MAX_LEN = rand() % (SUM_LEN - N50 + 1) + N50;
    int TMP_SUM = MAX_LEN;
    contig_list[0] = MAX_LEN;
    int count = 1;

    while (contig_list[count - 1] > N50 && count < TOT_SEQS) {
        int next_contig = rand() % contig_list[count - 1] + 1;
        contig_list[count++] = next_contig;
        TMP_SUM += next_contig;
    }

    if (count < TOT_SEQS && contig_list[count - 1] > N50) {
        contig_list[count - 1] = N50;
    }

    while (count < TOT_SEQS) {
        int remaining_sum = SUM_LEN - TMP_SUM;
        if (remaining_sum <= 0) break;
        int next_contig = rand() % MAX_LEN + 1;
        next_contig = (remaining_sum < next_contig) ? remaining_sum : next_contig;
        contig_list[count++] = next_contig;
        TMP_SUM += next_contig;
    }

    return contig_list;
}


void generate_sequence(int length, char *sequence) {
    memset(sequence, 'A', length);
    sequence[length] = '\0';
}

void write_fasta(char **sequences, int *lengths, int num_seqs, const char *outfile) {
    FILE *f = fopen(outfile, "w");
    if (!f) {
        fprintf(stderr, "Unable to write to file: %s\n", outfile);
        return;
    }

    for (int i = 0; i < num_seqs; i++) {
        fprintf(f, ">seq%d\n", i + 1);
        for (int j = 0; j < lengths[i]; j += 60) {
            int line_length = (lengths[i] - j < 60) ? lengths[i] - j : 60;
            fprintf(f, "%.*s\n", line_length, sequences[i] + j);
        }
    }

    fclose(f);
}

void write_fastq(char **sequences, int *lengths, int num_seqs, const char *outfile) {
    FILE *f = fopen(outfile, "w");
    if (!f) {
        fprintf(stderr, "Unable to write to file: %s\n", outfile);
        return;
    }

    for (int i = 0; i < num_seqs; i++) {
        fprintf(f, "@seq%d\n%s\n+\n", i + 1, sequences[i]);
        for (int j = 0; j < lengths[i]; j++) {
            fputc('I', f);  // Dummy quality score
        }
        fputc('\n', f);
    }

    fclose(f);
}

void free_resources(char **sequences, int *lengths, int num_seqs) {
    for (int i = 0; i < num_seqs; i++) {
        free(sequences[i]);
    }
    free(sequences);
    free(lengths);
}