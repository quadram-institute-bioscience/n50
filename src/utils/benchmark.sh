#!/usr/bin/env bash
set -euo pipefail

# Get input directory
INPUT_DIR="${1:-benchmark-data}"

if [[ ! -e "$INPUT_DIR" ]]; then
  echo "Input directory not found: $INPUT_DIR"
  exit 1
fi

OUT_DIR="$INPUT_DIR"/benchmark
FILE_DIR="$INPUT_DIR"/files
mkdir -p "$OUT_DIR"
mkdir -p "$FILE_DIR"


for FILE in "${INPUT_DIR}"/*.{fasta,fastq,fasta.gz,fastq.gz};
do
	BASE=$(basename $FILE  | sed 's/.fast.//' | sed 's/.gz$//')
        EXT=$(basename "$FILE" | sed 's/.*\.\(fasta\|fastq\(\.gz\)\?\)$/\1/')
	echo "$BASE | $EXT"
	if [[ ! -e "$FILE_DIR"/seqkit_${BASE}_${EXT}_output.txt ]]; then
          echo " - evaluating statistics"
	  n50 --format tsv "$FILE" > "$FILE_DIR"/n50v1_${BASE}_${EXT}_output.txt
	  if [[ $EXT == ".fasta" ]]; then
	    assembly_stats "$FILE" > "$FILE_DIR"/asstats_${BASE}_${EXT}_output.txt
          fi
	  bin/n50 "$FILE"  > "$FILE_DIR"/n50v2_${BASE}_${EXT}_output.txt
	  seqkit stats -T --all "$FILE" > "$FILE_DIR"/seqkit_${BASE}_${EXT}_output.txt
	fi
done

for FILE in "${INPUT_DIR}"/*.{fasta,fastq,fasta.gz,fastq.gz};
do
	BASE=$(basename $FILE  | sed 's/.fast.//' | sed 's/.gz$//')
        EXT=$(basename "$FILE" | sed 's/.*\.\(fasta\|fastq\(\.gz\)\?\)$/\1/')
	echo "$BASE | $EXT"
	if [[ $EXT == ".fasta" ]]; then
	  hyperfine --warmup 1 --max-runs 5 --export-markdown "$OUT_DIR"/${EXT}_${BASE}.md --export-csv "$OUT_DIR"/${EXT}_${BASE}.csv \
		-n seqkit "seqkit stats -T --all $FILE" \
		-n seqfu  "seqfu stats $FILE" \
		-n astats "assembly_stats $FILE" \
		-n n50_v2 "bin/n50 $FILE"
	else
	  hyperfine --warmup 1 --max-runs 5 --export-markdown "$OUT_DIR"/${EXT}_${BASE}.md --export-csv "$OUT_DIR"/${EXT}_${BASE}.csv \
		-n seqkit "seqkit stats -T --all $FILE" \
		-n seqfu  "seqfu stats $FILE" \
		-n n50_v2 "bin/n50 $FILE"
	fi
done
