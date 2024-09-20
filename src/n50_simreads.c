#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#define MAX_ARGS 100 // Maximum number of arguments
#define MAX_PATH 1024 // Maximum path length
/*
    n50 suite - simulate reads
    Usage: n50_simreads [--fasta|--fastq] -o OUTDIR ARGS

    ARGS format: COUNT*SIZE
    SIZE format: [0-9]+[KMG]?

    A program to simulate reads of given sizes and counts.
    Output is written to OUTDIR in FASTA or FASTQ format.
    Filename format: N50_TOTALSEQS_TOTLENGTH.fasta/fastq
*/
// Function to generate a random DNA sequence
void generate_sequence(char *seq, int length) {
    const char bases[] = "ACGTactAC"; // 70% AT, 30% CG
    for (int i = 0; i < length; i++) {
        seq[i] = bases[rand() % 9];
    }
    seq[length] = '\0';
}

// Function to generate a random quality string
void generate_quality(char *qual, int length) {
    for (int i = 0; i < length; i++) {
        qual[i] = 33 + (rand() % 41); // ASCII 33 to 73
    }
    qual[length] = '\0';
}

// Function to parse size string (e.g., "1kb", "2Mb")
int parse_size(const char *size_str) {
    int size = atoi(size_str);
    char suffix = toupper(size_str[strlen(size_str) - 1]);
    
    switch (suffix) {
        case 'K': size *= 1000; break;
        case 'M': size *= 1000000; break;
        case 'G': size *= 1000000000; break;
    }
    
    return size;
}

// Comparison function for qsort
int compare_ints(const void *a, const void *b) {
    return (*(int*)b - *(int*)a);
}

// Function to calculate N50
int calculate_n50(const int *lengths, int num_seqs, long long *total_length) {
    int *sorted_lengths = malloc(sizeof(int) * num_seqs);
    memcpy(sorted_lengths, lengths, sizeof(int) * num_seqs);
    
    qsort(sorted_lengths, num_seqs, sizeof(int), compare_ints);
    
    *total_length = 0;
    for (int i = 0; i < num_seqs; i++) {
        *total_length += sorted_lengths[i];
    }
    
    long long cumulative_length = 0;
    int n50 = -1;
    for (int i = 0; i < num_seqs; i++) {
        cumulative_length += sorted_lengths[i];
        if (cumulative_length >= *total_length / 2) {
            n50 = sorted_lengths[i];
            break;
        }
    }
    
    free(sorted_lengths);
    return n50;
}

int main(int argc, char *argv[]) {
    if (argc < 5) {
        fprintf(stderr, "Usage: %s [--fasta|--fastq] -o OUTDIR ARGS\n", argv[0]);
        return 1;
    }

    int is_fastq = 0;
    const char *format = is_fastq ? "FASTQ" : "FASTA";
    char *outdir = NULL;

    for (int i = 1; i < argc - 1; i++) {
        if (strcmp(argv[i], "--fastq") == 0) {
            is_fastq = 1;
        } else if (strcmp(argv[i], "--fasta") == 0) {
            is_fastq = 0;
        } else if (strcmp(argv[i], "-o") == 0) {
            outdir = argv[++i];
        }
    }

    if (!outdir) {
        fprintf(stderr, "Output directory not specified. Use -o OUTDIR.\n");
        return 1;
    }

    // Create output directory if it doesn't exist
    struct stat st = {0};
    if (stat(outdir, &st) == -1) {
        mkdir(outdir, 0700);
    }

    srand(1);

    int total_seqs = 0;
    int *lengths = NULL;
    int lengths_capacity = 0;    // Capacity of lengths array

    for (int i = 4; i < argc; i++) {
        char *arg = argv[i];
        char *count_str = strtok(arg, "*");
        char *size_str = strtok(NULL, "*");

        if (!count_str || !size_str) {
            fprintf(stderr, "Invalid argument format: %s\n", arg);
            continue;
        }

        int count = atoi(count_str);
        int size = parse_size(size_str);
        // print to stderr the count and size
        fprintf(stderr, "# Write seqs: %d, of size: %d\n", count, size);
        // Reallocate lengths array if necessary
        if (total_seqs + count > lengths_capacity) {
            lengths_capacity = total_seqs + count;
            lengths = realloc(lengths, lengths_capacity * sizeof(int));
            if (!lengths) {
                fprintf(stderr, "Memory allocation failed.\n");
                return 1;
            }
        }

        for (int j = 0; j < count; j++) {
            lengths[total_seqs++] = size;
            fprintf(stderr, " + Seq %d: %d\n", total_seqs, size);
        }
    }

    long long total_length;
    int n50 = calculate_n50(lengths, total_seqs, &total_length);

    // print N50, total seqs and total length to STDERR
    fprintf(stderr, "Format:\t%s\nN50:\t%d\nTot seqs:\t%d\nTot len:\t%lld\n", format, n50, total_seqs, total_length);
    char filename[MAX_PATH];
    snprintf(filename, MAX_PATH, "%s/%d_%d_%lld.%s", outdir, n50, total_seqs, total_length, is_fastq ? "fastq" : "fasta");

    FILE *outfile = fopen(filename, "w");
    if (!outfile) {
        fprintf(stderr, "Failed to open output file: %s\n", filename);
        free(lengths);
        return 1;
    }

    char *sequence = malloc(n50 + 1);
    char *quality = NULL;
    if (is_fastq) {
        quality = malloc(n50 + 1);
    }

    for (int i = 0; i < total_seqs; i++) {
        generate_sequence(sequence, lengths[i]);
        
        if (is_fastq) {
            generate_quality(quality, lengths[i]);
            fprintf(outfile, "@Simulated_read_%d len=%d\n%s\n+\n%s\n", i+1, lengths[i], sequence, quality);
        } else {
            fprintf(outfile, ">Simulated_read_%d len=%d\n%s\n", i+1, lengths[i],  sequence);
        }
    }

    fclose(outfile);
    free(sequence);
    if (quality) free(quality);
    free(lengths);

    printf("Output written to: %s\n", filename);

    return 0;
}