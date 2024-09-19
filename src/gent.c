#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <pthread.h>

#define MAX_SEQ_NAME_LEN 256
#define MAX_FILENAME_LEN 256
#define MAX_THREADS 16

// Thread-related structures
typedef struct {
    char **sequences;
    int *lengths;
    int start_index;
    int end_index;
} ThreadData;


// Function prototypes
long long calculate_n50(const int *lengths, int num_seqs, long long *total_length);
void *generate_sequences(void *arg);
void gen_ctg_len(int min_len, int max_len, int num_seqs, int *lengths);
void generate_sequence(int length, char *sequence);
void free_resources(char **sequences, int *lengths, int num_seqs);
void write_sequences(char **sequences, int *lengths, int num_seqs, const char *outfile, const char *format);

int main(int argc, char *argv[]) {
    if (argc != 9) {
        fprintf(stderr, "Usage: %s <min_seqs> <max_seqs> <min_len> <max_len> <tot_files> <format> <outdir> <num_threads>\n", argv[0]);
        return 1;
    }

    int min_seqs = atoi(argv[1]);
    int max_seqs = atoi(argv[2]);
    int min_len = atoi(argv[3]);
    int max_len = atoi(argv[4]);
    int tot_files = atoi(argv[5]);
    char *format = argv[6];
    char *outdir = argv[7];
    int num_threads = atoi(argv[8]);

    // Input validation
    if (min_seqs <= 0 || max_seqs <= 0 || min_len <= 0 || max_len <= 0 || tot_files <= 0 || num_threads <= 0 || num_threads > MAX_THREADS) {
        fprintf(stderr, "Error: Invalid input parameters.\n");
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
    fprintf(stderr, "Parameters: min_seqs=%d, max_seqs=%d, min_len=%d, max_len=%d, tot_files=%d, format=%s, outdir=%s, threads=%d\n",
            min_seqs, max_seqs, min_len, max_len, tot_files, format, outdir, num_threads);

    int seed = 42;  // Fixed seed for reproducibility
    srand(seed);

    for (int i = 1; i <= tot_files; i++) {
        int num_seqs = rand() % (max_seqs - min_seqs + 1) + min_seqs;
        int *contig_lengths = malloc(num_seqs * sizeof(int));
        if (!contig_lengths) {
            fprintf(stderr, "Memory allocation failed for contig_lengths\n");
            return 1;
        }
        fprintf(stderr, "%d/%d %s file (%d num seqs)", i, tot_files, format, num_seqs);
        gen_ctg_len(min_len, max_len, num_seqs, contig_lengths);
        fprintf(stderr, ", writing:\n");
        long long total_length;
        int N50 = calculate_n50(contig_lengths, num_seqs, &total_length);
        fprintf(stderr, "\tTotal length: %lld\n", total_length);
        fprintf(stderr, "\tN50: %d\n", N50);
        char outfile[MAX_FILENAME_LEN];
        snprintf(outfile, sizeof(outfile), "%s/%d_%d_%lld.%s", outdir, N50, num_seqs, total_length, format);

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
        }

        // Set up threading for sequence generation
        pthread_t threads[num_threads];
        ThreadData thread_data[num_threads];

        int seqs_per_thread = num_seqs / num_threads;
        int remaining_seqs = num_seqs % num_threads;

        for (int t = 0; t < num_threads; t++) {
            thread_data[t].sequences = sequences;
            thread_data[t].lengths = contig_lengths;
            thread_data[t].start_index = t * seqs_per_thread + (t < remaining_seqs ? t : remaining_seqs);
            thread_data[t].end_index = (t + 1) * seqs_per_thread + (t < remaining_seqs ? t + 1 : remaining_seqs);

            if (pthread_create(&threads[t], NULL, generate_sequences, (void *)&thread_data[t]) != 0) {
                fprintf(stderr, "Failed to create thread %d\n", t);
                free_resources(sequences, contig_lengths, num_seqs);
                return 1;
            }
        }

        // Wait for all threads to complete
        for (int t = 0; t < num_threads; t++) {
            pthread_join(threads[t], NULL);
        }

        // Write sequences to file (single-threaded)
        write_sequences(sequences, contig_lengths, num_seqs, outfile, format);

        fprintf(stderr, "  [Done: %s]\n", outfile);
        free_resources(sequences, contig_lengths, num_seqs);
    }

    return 0;
}


void gen_ctg_len(int min_len, int max_len, int num_seqs, int *lengths) {
    // calculate the even distance between num_seqs sequences spanning from min_len to max_len
    long long STEP = (max_len - min_len) / num_seqs;
    for (int i = 0; i < num_seqs; i++) {
        lengths[i] = STEP * i + min_len;
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
    
    // calculate int STEP defined as MAX_LEN
    int MAX_LEN = rand() % (SUM_LEN - N50 + 1) + N50;
    int TMP_SUM = MAX_LEN;
    contig_list[0] = MAX_LEN;
    int count = 1;

    while (contig_list[count - 1] > N50 && count < TOT_SEQS) {
        int next_contig = rand() % (int)(contig_list[count - 1] * 0.9) + 1;
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

void *generate_sequences(void *arg) {
    ThreadData *data = (ThreadData *)arg;
    for (int i = data->start_index; i < data->end_index; i++) {
        generate_sequence(data->lengths[i], data->sequences[i]);
    }
    return NULL;
}

void write_sequences(char **sequences, int *lengths, int num_seqs, const char *outfile, const char *format) {
    FILE *f = fopen(outfile, "w");
    if (!f) {
        fprintf(stderr, "Unable to write to file: %s\n", outfile);
        return;
    }

    if (strcmp(format, "fasta") == 0) {
        for (int i = 0; i < num_seqs; i++) {
            fprintf(f, ">seq_%d\n", i + 1);
            for (int j = 0; j < lengths[i]; j += 60) {
                int line_length = (lengths[i] - j < 60) ? lengths[i] - j : 60;
                fprintf(f, "%.*s\n", line_length, sequences[i] + j);
            }
        }
    } else {
        for (int i = 0; i < num_seqs; i++) {
            fprintf(f, "@seq%d\n%s\n+\n", i + 1, sequences[i]);
            for (int j = 0; j < lengths[i]; j++) {
                fputc('I', f);  // Dummy quality score
            }
            fputc('\n', f);
        }
    }

    fclose(f);
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