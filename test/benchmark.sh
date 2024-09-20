#!/usr/bin/env bash
set -euo pipefail
SELF_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PARENT_DIR="$(dirname "$SELF_DIR")"


# Get argument: INPUT_DIR
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 INPUT_DIR"
    exit 1
fi
INPUT_DIR=$1

# Check binary presence
for BIN in bin/n50 seqfu seqkit hyperfine;
do
  # Check binary presence with command -v
  if ! command -v "$BIN" &> /dev/null; then
    echo "Binary not found: $BIN"
    exit 1
  else
    echo "OK: Binary found: $BIN"
  fi
done

for FILE in "$INPUT_DIR"/*;
do 
 echo "== TEST $FILE"
 BASE=$(basename "$FILE" | sed 's/\./_/g')
 hyperfine --warmup 1 --max-runs 5 \
  --export-csv "$PARENT_DIR"/test/benchmark/single_"$BASE".csv \
  --export-markdown  "$PARENT_DIR"/test/benchmark/single_"$BASE".md \
   -n "n50" -n "seqfu" -n "seqkit" \
   "bin/n50 $FILE" \
   "seqfu stats $FILE" \
   "seqkit stats --all $FILE"
done