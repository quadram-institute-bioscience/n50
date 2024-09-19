#!/bin/bash
PARENT_DIR=$(dirname $(dirname $(readlink -f $0)))

mkdir -p $PARENT_DIR/test/benchmark
mkdir -p $PARENT_DIR/test/local
mkdir -p $PARENT_DIR/test/sim

for bin in seqfu seqkit ./bin/n50 ./bin/n50o; do
    if ! command -v $bin &> /dev/null; then
        echo "Error: $bin could not be found."
        echo "Please make sure that the binary is compiled and in the bin/ directory."
        exit 1
    fi
done

set -euo pipefail
# Let's if the subdirectories bin/  and test/ exist
if [ ! -d "$PARENT_DIR/bin" ] || [ ! -d "$PARENT_DIR/test" ] || [ ! -d "$PARENT_DIR/test/local" ]; then
    echo "Error: $PARENT_DIR/bin/ or $PARENT_DIR/test/ not found."
    echo "Please run the script from the root directory of the project."
    exit 1
fi


# the test/local directory contains files of different sizes and extensions
# *.fasta *.fastq *.fasta.gz and *.fastq.gz

ALL_FASTA=$(find $PARENT_DIR/test/local -type f -name "*.fasta" | xargs echo)
ALL_FASTQ=$(find $PARENT_DIR/test/local -type f -name "*.fastq" | xargs echo)
ALL_FASTA_GZ=$(find $PARENT_DIR/test/local -type f -name "*.fasta.gz" | xargs echo)
ALL_FASTQ_GZ=$(find $PARENT_DIR/test/local -type f -name "*.fastq.gz" | xargs echo)

# Files < 50Mb
FILES_SMALL=$(find $PARENT_DIR/test/local -type f -name "*.fast?" -size -50M | xargs echo)
# File > 500Mb
FILES_BIG=$(find $PARENT_DIR/test/local -type f -name "*.fast?" -size +500M | xargs echo)

# Count how many files we have
echo "FASTA files: $(echo $ALL_FASTA | wc -w)"
echo "FASTQ files: $(echo $ALL_FASTQ | wc -w)"
echo "FASTA GZ files: $(echo $ALL_FASTA_GZ | wc -w)"
echo "FASTQ GZ files: $(echo $ALL_FASTQ_GZ | wc -w)"
hyperfine --warmup 1 --max-runs 5 \
    --export-csv $PARENT_DIR/test/benchmark/ALL_FASTA.csv \
    "$PARENT_DIR/bin/n50 $ALL_FASTA" \
    "seqfu stats $ALL_FASTA" \
    "seqkit stats $ALL_FASTA" \
    "$PARENT_DIR/bin/n50o $ALL_FASTA"

hyperfine --warmup 1 --max-runs 5 \
    --export-csv $PARENT_DIR/test/benchmark/ALL_FASTA_GZ.csv \
    "$PARENT_DIR/bin/n50 $ALL_FASTA_GZ" \
    "seqfu stats $ALL_FASTA_GZ" \
    "seqkit stats $ALL_FASTA_GZ" \
    "$PARENT_DIR/bin/n50o $ALL_FASTA_GZ"

hyperfine --warmup 1 --max-runs 5 \
    --export-csv $PARENT_DIR/test/benchmark/ALL_FASTQ.csv \
    "$PARENT_DIR/bin/n50 $ALL_FASTQ" \
    "seqfu stats $ALL_FASTQ" \
    "seqkit stats $ALL_FASTQ" \
    "$PARENT_DIR/bin/n50o $ALL_FASTQ"

hyperfine --warmup 1 --max-runs 5 \
    --export-csv $PARENT_DIR/test/benchmark/ALL_FASTQ_GZ.csv \
    "$PARENT_DIR/bin/n50 $ALL_FASTQ_GZ" \
    "seqfu stats $ALL_FASTQ_GZ" \
    "seqkit stats $ALL_FASTQ_GZ" \
    "$PARENT_DIR/bin/n50o $ALL_FASTQ_GZ"


hyperfine --warmup 1 --max-runs 5 \
    --export-csv $PARENT_DIR/test/benchmark/ALL_SMALL.csv \
    "$PARENT_DIR/bin/n50 $FILES_SMALL" \
    "seqfu stats $FILES_SMALL" \
    "seqkit stats $FILES_SMALL" \
    "$PARENT_DIR/bin/n50o $FILES_SMALL"

hyperfine --warmup 1 --max-runs 5 \
    --export-csv $PARENT_DIR/test/benchmark/ALL_BIG.csv \
    "$PARENT_DIR/bin/n50 $FILES_BIG" \
    "seqfu stats $FILES_BIG" \
    "seqkit stats $FILES_BIG" \
    "$PARENT_DIR/bin/n50o $FILES_BIG"