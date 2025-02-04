#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <limits.h>

#define MAX_ARGS 100 // Maximum number of arguments
#define MAX_PATH 1024 // Maximum path length

#define MAX_NUM_LENGTH 30  // Enough for 64-bit integers

// make const char bases[] = "ACGTactAC"; // 70% AT, 30% CG as global constant
const char bases[] = "ACGTactAC"; 
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
    if (!seq || length < 0) {
        return;
    }
    
    for (int i = 0; i < length; i++) {
        seq[i] = bases[rand() % (sizeof(bases) - 1)];  // -1 to exclude null terminator
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


char* num_to_str(long long number, char* out_str, size_t size) {
    if (!out_str || size < 2) return NULL; // Basic validation
    
    int is_negative = 0;
    int len = 0;
    long long abs_number;
    char temp[MAX_NUM_LENGTH]; // Temporary buffer for building the string

    // Handle negative numbers
    if (number < 0) {
        is_negative = 1;
        abs_number = -number;
    } else {
        abs_number = number;
    }

    // Convert number to string (reverse order)
    do {
        temp[len++] = abs_number % 10 + '0';
        abs_number /= 10;
    } while (abs_number > 0 && len < (int)size - 2); // Leave room for sign and null

    // Add commas
    for (int i = 3; i < len; i += 4) {
        if (len + 1 >= (int)size - 1) break;  // Prevent buffer overflow
        memmove(&temp[i + 1], &temp[i], len - i + 1);
        temp[i] = ',';
        len++;
    }

    // Add minus sign if negative
    if (is_negative) {
        if (len + 1 >= (int)size) len--;  // Make room if necessary
        memmove(&temp[1], &temp[0], len + 1);
        temp[0] = '-';
        len++;
    }

    // Reverse the string into the output buffer
    for (int i = 0; i < len; i++) {
        out_str[i] = temp[len - 1 - i];
    }
    out_str[len] = '\0';

    return out_str;
}

// Function to parse size string (e.g., "1kb", "2Mb")
// Safer parse_size with overflow checking
long long parse_size(const char *size_str) {
    char *endptr;
    long long size = strtoll(size_str, &endptr, 10);
    if (size < 0) {
        fprintf(stderr, "Error: Negative size not allowed\n");
        return -1;
    }
    
    if (endptr == size_str) {
        fprintf(stderr, "Error: Invalid number format\n");
        return -1;
    }
    
    char suffix = toupper(*endptr);
    long long multiplier = 1;
    
    switch (suffix) {
        case 'K': multiplier = 1000LL; break;
        case 'M': multiplier = 1000000LL; break;
        case 'G': multiplier = 1000000000LL; break;
        case '\0': return size; // No suffix
        default:
            fprintf(stderr, "Error: Invalid suffix '%c'\n", suffix);
            return -1;
    }
    
    // Check for overflow
    if (size > LLONG_MAX / multiplier) {
        fprintf(stderr, "Error: Size too large\n");
        return -1;
    }
    
    size *= multiplier;
    return size;
}

// Comparison function for qsort
int compare_ints(const void *a, const void *b) {
    return (*(int*)b - *(int*)a);
}
int compare_longs(const void *a, const void *b) {
    long long va = *(const long long*)a;
    long long vb = *(const long long*)b;
    return (va > vb) - (va < vb);
}
// Function to calculate N50
long long calculate_n50(const long long *lengths, long long num_seqs, long long *total_length) {
    long long *sorted_lengths = malloc(sizeof(long long) * num_seqs);
    memcpy(sorted_lengths, lengths, sizeof(long long) * num_seqs);
    
    qsort(sorted_lengths, num_seqs, sizeof(long long), compare_longs);
    
    *total_length = 0;
    for (long long i = 0; i < num_seqs; i++) {
        *total_length += sorted_lengths[i];
    }
    
    long long cumulative_length = 0;
    long long n50 = -1;
    for (long long i = 0; i < num_seqs; i++) {
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
        fprintf(stderr, "Usage: %s [--fasta|--fastq] -o OUTDIR [-p PREFIX] ARGS\n", argv[0]);
        fprintf(stderr, "ARGS format: COUNT*SIZE\n");
        return 1;
    } else if (argc > MAX_ARGS) {
        fprintf(stderr, "Too many arguments. Maximum is %d.\n", MAX_ARGS);
        return 1;
    }
    long long max_seq_length = 0;
    int is_fastq = 0;
    int verbose = 0;
    char *outdir = NULL;
    char *prefix = NULL;
    const char *format = NULL;
    char *sequence = NULL;
    char *quality = NULL;
    long long *lengths = NULL;
    FILE *outfile = NULL;
    int ret = 1;  // Default to error
    long long total_seqs = 0;
    long long lengths_capacity = 0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--fastq") == 0) {
            is_fastq = 1;
        } else if (strcmp(argv[i], "--fasta") == 0) {
            is_fastq = 0;
        } else if (strcmp(argv[i], "-o") == 0) {
            outdir = argv[++i];
        } else if (strcmp(argv[i], "--verbose") == 0) {
            verbose = 1;
            fprintf(stderr, "Verbose mode enabled.\n");
        } else if (strcmp(argv[i], "-p") == 0) {
            prefix = argv[++i];
        } else if (!strchr(argv[i], '*')) {
            fprintf(stderr, "Invalid argument: %s\n", argv[i]);
            goto cleanup;
        }
    }
    
    format = is_fastq ? "FASTQ" : "FASTA";

    if (!outdir) {
        fprintf(stderr, "Output directory not specified. Use -o OUTDIR.\n");
        goto cleanup;
    }
    if (!prefix) {
        prefix = "";
    }

    // Create output directory if it doesn't exist
    struct stat st = {0};
    if (stat(outdir, &st) == -1) {
        mkdir(outdir, 0700);
    }

    srand(1);

    for (int i = 0; i < argc; i++) {
        if (!strchr(argv[i], '*')) {
            continue;
        }
        char *arg = strdup(argv[i]);  // Create a copy since strtok modifies the string
        if (!arg) {
            fprintf(stderr, "Memory allocation failed.\n");
            goto cleanup;
        }
        
        char *count_str = strtok(arg, "*");
        char *size_str = strtok(NULL, "*");

        if (!count_str || !size_str) {
            fprintf(stderr, "Invalid argument format: %s\n", arg);
            free(arg);
            continue;
        }

        long long count = atoll(count_str);
        long long size = parse_size(size_str);
        free(arg);


    if (size > max_seq_length) {
        max_seq_length = size;
    }

    if (total_seqs + count > lengths_capacity) {
        lengths_capacity = total_seqs + count;
        long long *new_lengths = realloc(lengths, lengths_capacity * sizeof(long long));
        if (!new_lengths) {
            fprintf(stderr, "Memory allocation failed.\n");
            goto cleanup;
        }
        lengths = new_lengths;
    }

    for (long long j = 0; j < count; j++) {
        lengths[total_seqs++] = size;
    }
    }

    long long total_length;
    long long n50 = calculate_n50(lengths, total_seqs, &total_length);

    char n50_str[MAX_NUM_LENGTH];
    char seqs_str[MAX_NUM_LENGTH];
    char len_str[MAX_NUM_LENGTH];

    fprintf(stderr, "\n------\nMode:\t%s\nPrefix:\t%s\nFormat:\t%s\nN50:\t%s\nTot seqs:\t%s\nTot len:\t%s\n------\n",
        verbose ? "verbose" : "standard", prefix, format,
        num_to_str(n50, n50_str, sizeof(n50_str)),
        num_to_str(total_seqs, seqs_str, sizeof(seqs_str)),
        num_to_str(total_length, len_str, sizeof(len_str)));

    char filename[MAX_PATH];
    int path_len = snprintf(filename, sizeof(filename), 
                          "%s/%s%lld_%lld_%lld.%s", 
                          outdir, prefix, n50, total_seqs, total_length,
                          is_fastq ? "fastq" : "fasta");
    
    if ((size_t)path_len >= sizeof(filename)) {
        fprintf(stderr, "Output path too long\n");
        goto cleanup;
    }

    outfile = fopen(filename, "w");
    if (!outfile) {
        fprintf(stderr, "Failed to open output file: %s\n", filename);
        goto cleanup;
    }

    sequence = malloc(max_seq_length + 1);  // +1 for null terminator
    if (!sequence) {
        fprintf(stderr, "Failed to allocate sequence buffer\n");
        goto cleanup;
    }


    if (is_fastq) {
        quality = malloc(max_seq_length + 1);  // +1 for null terminator
        if (!quality) {
            fprintf(stderr, "Failed to allocate quality buffer\n");
            goto cleanup;
        }
    }

    int last_length = 0;
    for (int i = 0; i < total_seqs; i++) {
        if (lengths[i] != last_length) {
            last_length = lengths[i];
        }

        generate_sequence(sequence, lengths[i]);
        if (verbose && i % 1000 == 0) {
            fprintf(stderr, " Generating seq #%d (%lld bp)\r", i, lengths[i]);
        }
        
        if (is_fastq) {
            generate_quality(quality, lengths[i]);
            fprintf(outfile, "@Simulated_read_%d len=%lld\n%s\n+\n%s\n", i+1, lengths[i], sequence, quality);
        } else {
            fprintf(outfile, ">Simulated_read_%d len=%lld\n%s\n", i+1, lengths[i], sequence);
        }
    }

    fprintf(stderr, "\n");
    printf("Output written to: %s\n", filename);
    ret = 0;  // Success

cleanup:
    if (outfile) fclose(outfile);
    free(sequence);
    free(quality);
    free(lengths);
    return ret;
}