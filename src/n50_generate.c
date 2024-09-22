#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <unistd.h>
#include <libgen.h>

#define MAX_LINE_LENGTH 1024
#define MAX_COMMAND_LENGTH 4096

void print_usage(const char *program_name) {
    fprintf(stderr, "Usage: %s -i INPUTFILE -o OUTDIR [-f FORMAT] [-s PATH]\n", program_name);
    fprintf(stderr, "  -i INPUTFILE   Input file path\n");
    fprintf(stderr, "  -o OUTDIR      Output directory\n");
    fprintf(stderr, "  -f FORMAT      Optional: Output format (FASTQ by default, FASTA also supported)\n");
    fprintf(stderr, "  -s PATH        Optional: Path to n50_simreads executable\n");
}

char *get_executable_path(const char *argv0) {
    char *path = strdup(argv0);
    if (!path) {
        perror("Memory allocation failed");
        exit(EXIT_FAILURE);
    }
    char *dir = dirname(path);
    char *full_path = malloc(strlen(dir) + strlen("/n50_simreads") + 1);
    if (!full_path) {
        perror("Memory allocation failed");
        free(path);
        exit(EXIT_FAILURE);
    }
    sprintf(full_path, "%s/n50_simreads", dir);
    free(path);
    return full_path;
}

void process_input_file(const char *input_file, char *generated_string, size_t max_length) {
    FILE *file = fopen(input_file, "r");
    if (!file) {
        perror("Error opening input file");
        exit(EXIT_FAILURE);
    }

    char line[MAX_LINE_LENGTH];
    int line_number = 0;
    size_t current_length = 0;

    while (fgets(line, sizeof(line), file)) {
        line_number++;
        if (line_number == 1) continue; // Skip header

        char *token = strtok(line, ",");
        if (!token) continue;

        int length = atoi(token);
        token = strtok(NULL, ",");
        if (!token) continue;

        int count = atoi(token);

        if (length <= 0 || count <= 0) {
            fprintf(stderr, "Warning: Invalid data on line %d, skipping\n", line_number);
            continue;
        }

        int written = snprintf(generated_string + current_length, max_length - current_length,
                               "%d*%d ", count, length);

        if (written < 0 || written >= max_length - current_length) {
            fprintf(stderr, "Error: Generated string too long\n");
            fclose(file);
            exit(EXIT_FAILURE);
        }

        current_length += written;
    }

    if (current_length > 0) {
        generated_string[current_length - 1] = '\0'; // Remove trailing space
    }

    fclose(file);
}

void run_n50_simreads(const char *format, const char *outdir, const char *generated_string, const char *n50_simreads_path) {
    char command[MAX_COMMAND_LENGTH];
    snprintf(command, sizeof(command), "%s --%s -o %s %s",
             n50_simreads_path, format, outdir, generated_string);

    printf("Executing command: %s\n", command);

    int result = system(command);
    if (result == -1) {
        perror("Error executing n50_simreads command");
        exit(EXIT_FAILURE);
    } else if (WIFEXITED(result) && WEXITSTATUS(result) != 0) {
        fprintf(stderr, "n50_simreads command failed with exit status %d\n", WEXITSTATUS(result));
        exit(EXIT_FAILURE);
    }
}

void calculate_stats(const char *input_file) {
    FILE *file = fopen(input_file, "r");
    if (!file) {
        perror("Error opening input file for stats calculation");
        exit(EXIT_FAILURE);
    }

    char line[MAX_LINE_LENGTH];
    int line_number = 0;
    long long total_reads = 0;
    int max_length = 0;

    while (fgets(line, sizeof(line), file)) {
        line_number++;
        if (line_number == 1) continue; // Skip header

        char *token = strtok(line, ",");
        if (!token) continue;

        int length = atoi(token);
        token = strtok(NULL, ",");
        if (!token) continue;

        int count = atoi(token);

        if (length <= 0 || count <= 0) {
            fprintf(stderr, "Warning: Invalid data on line %d, skipping\n", line_number);
            continue;
        }

        total_reads += count;
        if (length > max_length) {
            max_length = length;
        }
    }

    fclose(file);

    printf("Total number of reads: %lld\n", total_reads);
    printf("Maximum read length: %d\n", max_length);
}

int main(int argc, char *argv[]) {
    char *input_file = NULL;
    char *outdir = NULL;
    char *format = "FASTQ";
    char *n50_simreads_path = NULL;

    int opt;
    while ((opt = getopt(argc, argv, "i:o:f:s:h")) != -1) {
        switch (opt) {
            case 'i':
                input_file = optarg;
                break;
            case 'o':
                outdir = optarg;
                break;
            case 'f':
                format = optarg;
                break;
            case 's':
                n50_simreads_path = optarg;
                break;
            case 'h':
                print_usage(argv[0]);
                exit(EXIT_SUCCESS);
            default:
                print_usage(argv[0]);
                exit(EXIT_FAILURE);
        }
    }

    if (!input_file || !outdir) {
        fprintf(stderr, "Error: Input file and output directory are required.\n");
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    // Convert format to uppercase
    for (char *p = format; *p; ++p) *p = toupper(*p);

    // Validate format
    if (strcmp(format, "FASTQ") != 0 && strcmp(format, "FASTA") != 0) {
        fprintf(stderr, "Error: Invalid format. Use FASTQ or FASTA.\n");
        exit(EXIT_FAILURE);
    }

    // Set default n50_simreads path if not provided
    if (!n50_simreads_path) {
        n50_simreads_path = get_executable_path(argv[0]);
    }

    // Process input file
    char generated_string[MAX_COMMAND_LENGTH];
    process_input_file(input_file, generated_string, sizeof(generated_string));

    // Run n50_simreads
    run_n50_simreads(format, outdir, generated_string, n50_simreads_path);

    // Calculate and print stats
    calculate_stats(input_file);

    // Free allocated memory
    if (n50_simreads_path != NULL && n50_simreads_path != get_executable_path(argv[0])) {
        free(n50_simreads_path);
    }

    return 0;
}