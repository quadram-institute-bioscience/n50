#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>
#include <unistd.h>

#define VERSION "1.9.3"

#define MAX_READS 1000000
#define BASES "ACGTactAC"
#define BASES_LEN 9
#define QUAL_MIN 33
#define QUAL_MAX 73

typedef struct {
    long count;
    long size;
} ReadSpec;

long parse_size(const char *str) {
    char *end;
    long size = strtol(str, &end, 10);
    if (*end) {
        switch (toupper(*end)) {
            case 'K': size *= 1000; break;
            case 'M': size *= 1000000; break;
            case 'G': size *= 1000000000; break;
            default:
                fprintf(stderr, "Invalid size suffix: %c\n", *end);
                return 0;
        }
    }
    return size;
}

int cmp_desc(const void *a, const void *b) {
    long sa = *(long*)a;
    long sb = *(long*)b;
    return sb - sa;
}

long compute_n50(long *sizes, int n) {
    if (n == 0) return 0; // Handle empty array case
    qsort(sizes, n, sizeof(long), cmp_desc);
    long total = 0;
    for (int i = 0; i < n; i++) total += sizes[i];
    long half = total / 2;
    long acc = 0;
    for (int i = 0; i < n; i++) {
        acc += sizes[i];
        if (acc >= half) return sizes[i];
    }
    return 0;
}

void rand_seq(char *buf, long len) {
    for (long i = 0; i < len; i++)
        buf[i] = BASES[rand() % BASES_LEN];
    buf[len] = '\0';
}

void rand_qual(char *buf, long len) {
    for (long i = 0; i < len; i++)
        buf[i] = QUAL_MIN + rand() % (QUAL_MAX - QUAL_MIN + 1);
    buf[len] = '\0';
}

void print_help(const char *prog) {    if (!prog) prog = "n50_simreads";    fprintf(stderr, "Usage: %s [--fasta|--fastq] -o OUTDIR [-p PREFIX] ARGS\n", prog);    fprintf(stderr, "ARGS format: COUNT*SIZE\n");    fprintf(stderr, "  --version  Show version number and exit\n");}

int main(int argc, char *argv[]) {
    int fasta = 1, verbose = 0;
    char *outdir = NULL, *prefix = "";
    ReadSpec *specs = malloc(MAX_READS * sizeof(ReadSpec));
    if (!specs) {
        fprintf(stderr, "Memory allocation failed for specs.\n");
        return 1;
    }
    int spec_count = 0;
    long total_reads = 0, total_bases = 0;

    srand(1);

    if (argc == 1) {
        print_help(argv[0]);
        return 0;
    }

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--fasta") == 0) fasta = 1;
        else if (strcmp(argv[i], "--fastq") == 0) fasta = 0;
        else if (strcmp(argv[i], "--verbose") == 0) verbose = 1;
        else if ((strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "--help") == 0)) {
            print_help(argv[0]);
            return 0;
        } else if (strcmp(argv[i], "--version") == 0) {
            printf("%s\n", VERSION);
            return 0;
        } else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) outdir = argv[++i];
        else if (strcmp(argv[i], "-p") == 0 && i + 1 < argc) prefix = argv[++i];
        else if (strchr(argv[i], '*')) {
            char *x = strdup(argv[i]);
            char *c = strchr(x, '*');
            if (c == NULL) { // Argument does not contain '*', treat as unknown
                fprintf(stderr, "Unknown option or argument: %s\n", argv[i]);
                free(x);
                return 1;
            }
            *c = '\0';
            long count = atol(x);
            long size = parse_size(c + 1);
            if (count <= 0 || size <= 0) {
                fprintf(stderr, "Invalid COUNT*SIZE: %s\n", argv[i]);
                free(x);
                free(specs);
                return 1;
            }
            specs[spec_count].count = count;
            specs[spec_count].size = size;
            spec_count++;
            total_reads += count;
            total_bases += count * size;
            free(x);
        } else {
            fprintf(stderr, "Unknown option or argument: %s\n", argv[i]);
            return 1;
        }
    }

    if (!outdir) {
        fprintf(stderr, "Output directory (-o) is required.\n");
        free(specs);
        return 1;
    }

    if (total_reads == 0) {
        fprintf(stderr, "Error: No read specifications provided. Please provide arguments in COUNT*SIZE format.\n");
        print_help(argv[0]);
        free(specs);
        return 1;
    }

    if (mkdir(outdir, 0755) != 0 && access(outdir, W_OK) != 0) {
        perror("Failed to create output directory");
        free(specs);
        return 1;
    }

    long *sizes = malloc(total_reads * sizeof(long));    if (!sizes) {        fprintf(stderr, "Memory allocation failed for sizes.\n");        free(specs);        return 1;    }

    long idx = 0;
    for (int i = 0; i < spec_count; i++)
        for (long j = 0; j < specs[i].count; j++)
            sizes[idx++] = specs[i].size;

    long n50 = compute_n50(sizes, total_reads);

    char filename[1024];
    snprintf(filename, sizeof(filename), "%s/%s%ld_%ld_%ld.%s", outdir, prefix, n50, total_reads, total_bases, fasta ? "fasta" : "fastq");

    FILE *out = fopen(filename, "w");
    if (!out) {
        perror("fopen");
        free(sizes);
        free(specs);
        return 1;
    }

    long max_size = 0;
    for (int i = 0; i < spec_count; i++) {
        if (specs[i].size > max_size) {
            max_size = specs[i].size;
        }
    }

    char *seq = malloc(max_size + 1);
    char *qual = malloc(max_size + 1);
    if (!seq || !qual) {
        fprintf(stderr, "Memory allocation failed for sequences.\n");
        fclose(out);
        free(sizes);
        free(specs);
        return 1;
    }

    idx = 0;
    for (int i = 0; i < spec_count; i++) {
        for (long j = 0; j < specs[i].count; j++) {
            rand_seq(seq, specs[i].size);
            if (fasta)
                fprintf(out, ">read%ld\n%s\n", idx++, seq);
            else {
                rand_qual(qual, specs[i].size);
                fprintf(out, "@read%ld\n%s\n+\n%s\n", idx++, seq, qual);
            }
            if (verbose && idx % 10000 == 0)
                fprintf(stderr, "Generated %ld reads...\n", idx);
        }
    }

    free(seq);
    free(qual);
    free(sizes);
    free(specs);
    fclose(out);

    if (verbose) {
        fprintf(stderr, "Output written to: %s\n", filename);
    }

    return 0;
}

