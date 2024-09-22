#!/bin/bash

# Check if required arguments are provided
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <input_file> <FORMAT> <OUTDIR>"
    exit 1
fi

INPUT_FILE=$1
FORMAT=$2
OUTDIR=$3

# Process the input file and generate the string
GENERATED_STRING=$(awk -F',' '
    NR>1 {
        if ($1 ~ /^[0-9]+$/ && $2 ~ /^[0-9]+$/ && $2 > 0) {
            printf "%d*%d ", $2, $1
        }
    }
' "$INPUT_FILE" | sed 's/ $//')

# Run the n50_simreads command with the generated string
n50_simreads --${FORMAT} -o ${OUTDIR} $GENERATED_STRING

# Print the command that was executed (for verification)
echo "Executed command: n50_simreads --${FORMAT} -o ${OUTDIR} $GENERATED_STRING"

# Calculate and print total number of reads
TOTAL_READS=$(awk -F',' 'NR>1 && $2 ~ /^[0-9]+$/ {sum += $2} END {print sum}' "$INPUT_FILE")
echo "Total number of reads: $TOTAL_READS"

# Find the maximum read length
MAX_LENGTH=$(awk -F',' 'NR>1 && $1 ~ /^[0-9]+$/ && $2 > 0 {max=$1} END {print max}' "$INPUT_FILE")
echo "Maximum read length: $MAX_LENGTH"
