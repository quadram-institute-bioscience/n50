#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#define MAX_LINE 1000000
#define NUM_BINS 16

// Function to determine the bin for a given length
int get_bin(int length) {
    int bins[] = {10, 100, 1000, 2500, 5000, 10000, 20000, 35000, 50000, 75000, 100000, 200000, 300000, 500000, 750000, 1000000};
    for (int i = 0; i < NUM_BINS; i++) {
        if (length <= bins[i]) {
            return i;
        }
    }
    return NUM_BINS - 1;  // For reads longer than the last bin
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <fastq_file>\n", argv[0]);
        return 1;
    }

    gzFile fp = gzopen(argv[1], "r");
    if (!fp) {
        fprintf(stderr, "Error: Could not open file %s\n", argv[1]);
        return 1;
    }

    char line[MAX_LINE];
    int counters[NUM_BINS] = {0};
    int line_count = 0;

    while (gzgets(fp, line, sizeof(line))) {
        line_count++;
        if (line_count % 4 == 2) {  // This is the sequence line
            int length = strlen(line) - 1;  // Subtract 1 to remove newline
            int bin = get_bin(length);
            counters[bin]++;
        }
    }

    gzclose(fp);

    // Print results
    printf("Bin,Number of Reads\n");
    int bins[] = {10, 100, 1000, 2500, 5000, 10000, 20000, 35000, 50000, 75000, 100000, 200000, 300000, 500000, 750000, 1000000};
    for (int i = 0; i < NUM_BINS; i++) {
        printf("%d,%d\n", bins[i], counters[i]);
    }

    return 0;
}