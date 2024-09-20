#!/usr/bin/env bash
set -euo pipefail
SELF_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PARENT_DIR="$(dirname "$SELF_DIR")"
BIN_DIR="$PARENT_DIR/bin"
# Check binary presence
for BIN in n50 n50_simreads gen;
do
    if [ ! -f "$BIN_DIR/$BIN" ]; then
        echo "Binary not found: $BIN"
        exit 1
    else
        echo "OK: Binary found: $BIN"
    fi
done

# get 1 argument to perform deep test
DEEP=0
if [ "$#" -eq 1 ]; then
    DEEP=1
fi

# Simulate reads
OUT_DIR="$PARENT_DIR/test/sim/"
if [ -d "$OUT_DIR" ]; then
    echo "Cleaning $OUT_DIR"
    # Remove all files if not empty already
    rm -f "$OUT_DIR"/* || true
fi
mkdir -p "$OUT_DIR"
COUNTER=0
# $COMPRESSOR = gzip or pigz if available
COMPRESSOR=$(which pigz)
if [ -z "$COMPRESSOR" ]; then
    COMPRESSOR=$(which gzip)
fi

for FORMAT in fasta fastq;
do
    COUNTER=$((COUNTER+4))
    # Small/Medium
    "${BIN_DIR}"/n50_simreads --${FORMAT} -o $OUT_DIR 100*250   20*1k   1*50k  2> /dev/null
    "${BIN_DIR}"/n50_simreads --${FORMAT} -o $OUT_DIR 100*1000  10*100k 1*10M  2> /dev/null
    # Large
    if [[ $FORMAT == "fastq" ]] && [[ $DEEP == 1 ]]; then
        "${BIN_DIR}"/n50_simreads --${FORMAT} -o $OUT_DIR 100*5k 1*1G  2> /dev/null
    fi
    ${COMPRESSOR} -k "$OUT_DIR"/*."$FORMAT"
done

I=0
for FILE in "$OUT_DIR"/*; 
do
    I=$((I+1))
    echo "$I/$COUNTER: Simulated reads: $FILE"
    OUTPUT=$("${BIN_DIR}"/n50 "$FILE" | sed 's/\t/,/g')
    # First field
    FILENAME=$(basename $(echo $OUTPUT | cut -d' ' -f1))
    EXP_N50=$(echo "${FILENAME}" | cut -f 1 -d '_' | cut -f 1 -d '.')
    EXP_SEQS=$(echo "${FILENAME}" | cut -f 2 -d '_' | cut -f 1 -d '.')
    EXP_BASES=$(echo "${FILENAME}" | cut -f 3 -d '_'  | cut -f 1 -d '.')
    REAL_N50=$(echo "$OUTPUT" | cut -d',' -f 5)
    REAL_SEQS=$(echo "$OUTPUT" | cut -d',' -f 4)
    REAL_BASES=$(echo "$OUTPUT" | cut -d',' -f 3)
    if [ "$EXP_N50" != "$REAL_N50" ]; then
        echo "ERROR: N50 mismatch: $EXP_N50 != $REAL_N50"
        exit 1
    else
        echo "OK: N50 match: $EXP_N50 == $REAL_N50"
    fi

    if [ "$EXP_SEQS" != "$REAL_SEQS" ]; then
        echo "ERROR: Sequences mismatch: $EXP_SEQS != $REAL_SEQS"
        exit 1
    else
        echo "OK: Sequences match: $EXP_SEQS == $REAL_SEQS"
    fi


    if [ "$EXP_BASES" != "$REAL_BASES" ]; then
        echo "ERROR: Bases mismatch: $EXP_BASES != $REAL_BASES"
        exit 1
    else
        echo "OK: Bases match: $EXP_BASES == $REAL_BASES"
    fi
    echo ""
done

# if hyperfine and seqkit and seqfu are available go on
if [ -z "$(which hyperfine)" ] || [ -z "$(which seqkit)" ] || [ -z "$(which seqfu)" ]; then
    echo "Skipping benchmark"
    exit 0
fi

mkdir -p "$PARENT_DIR"/test/benchmark/

hyperfine --warmup 1 --max-runs 9 \
  --export-csv "$PARENT_DIR"/test/benchmark/all_files.csv \
  -n "n50" -n "seqfu" -n "seqkit" \
  "${BIN_DIR}/n50 $OUT_DIR/*" \
  "seqfu stats $OUT_DIR/*" \
  "seqkit stats --all $OUT_DIR/*"
  
hyperfine --warmup 1 --max-runs 9 \
  --export-csv "$PARENT_DIR"/test/benchmark/all_fastq.csv \
    -n "n50" -n "seqfu" -n "seqkit" \
  "${BIN_DIR}/n50 $OUT_DIR/*.fastq" \
  "seqfu stats $OUT_DIR/*.fastq" \
  "seqkit stats --all $OUT_DIR/*.fastq"

hyperfine --warmup 1 --max-runs 9 \
  --export-csv "$PARENT_DIR"/test/benchmark/all_fasta.csv \
    -n "n50" -n "seqfu" -n "seqkit" \
  "${BIN_DIR}/n50 $OUT_DIR/*.fasta" \
  "seqfu stats $OUT_DIR/*.fasta" \
  "seqkit stats --all $OUT_DIR/*.fasta"

hyperfine --warmup 1 --max-runs 9 \
  --export-csv "$PARENT_DIR"/test/benchmark/all_fasta_gz.csv \
    -n "n50" -n "seqfu" -n "seqkit" \
  "${BIN_DIR}/n50 $OUT_DIR/*.fasta.gz" \
  "seqfu stats $OUT_DIR/*.fasta.gz" \
  "seqkit stats --all $OUT_DIR/*.fasta.gz"

hyperfine --warmup 1 --max-runs 9 \
  --export-csv "$PARENT_DIR"/test/benchmark/all_fastq_gz.csv \
    -n "n50" -n "seqfu" -n "seqkit" \
  "${BIN_DIR}/n50 $OUT_DIR/*.fastq.gz" \
  "seqfu stats $OUT_DIR/*.fastq.gz" \
  "seqkit stats --all $OUT_DIR/*.fastq.gz"

for FILE in "$OUT_DIR"/*;
do
    echo "Benchmarking single file: $FILE"
    BASE=$(basename $FILE | sed 's/\./_/g')
    hyperfine --warmup 1 --max-runs 9 \
      --export-csv $PARENT_DIR/test/benchmark/single_${BASE}.csv \
        -n "n50" -n "seqfu" -n "seqkit" \
        "${BIN_DIR}/n50 $FILE" \
        "seqfu stats $FILE" \
        "seqkit stats --all $FILE"
done