#!/bin/bash
set -euo pipefail

# Color definitions
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[1;34m'
NC='\033[0m' # No Color

header() {
    echo -e "${BLUE}==> $1${NC}"
}

success() {
    echo -e "${GREEN}✔ $1${NC}"
}

fail() {
    echo -e "${RED}✖ $1${NC}"
}

info() {
    echo -e "${YELLOW}→ $1${NC}"
}

# Start build
header "Building executables..."
make

# Version checks
V1=$(./bin/n50 --version)
V2=$(./bin/n50_simseqs --version)

header "Checking versions..."
if [[ $V1 =~ [0-9]+\.[0-9]+\.[0-9]+$ ]]; then
    success "Valid version: n50 $V1"
else
    fail "Invalid version: n50 $V1"
fi

if [[ $V2 =~ [0-9]+\.[0-9]+\.[0-9]+$ ]]; then
    success "Valid version: n50_simseqs $V2"
else
    fail "Invalid version: n50_simseqs $V2"
fi

# Run test
header "Running minimal N50 test..."
./bin/n50 --basename ./test/test.fa > ./test/output.tsv

expected_output="test.fa	3	34	18	12	4	1	52.94	11.33	4	18"
actual_output=$(tail -n 1 ./test/output.tsv)

if [ "$actual_output" == "$expected_output" ]; then
    success "N50 minimal test passed!"
else
    fail "N50 minimal test failed!"
    info "Expected: $expected_output"
    info "Actual:   $actual_output"
    exit 1
fi

# Simulate data
header "Generating synthetic sequences..."
OUTDIR="test-data"
mkdir -p "$OUTDIR"

for format in fasta fastq; do
    bin/n50_simseqs --${format} -o ${OUTDIR} -p test_ 1*20M 1*1M 10*9K 100*1K 1000*120
    bin/n50_simseqs --${format} -o ${OUTDIR} -p test_ 10*1M 100*2K 1000*120 2000*50
done

# Compress
header "Compressing outputs..."
for FILE in ${OUTDIR}/*.{fasta,fastq}; do
    if [[ ! -e ${FILE}.gz ]]; then
        info "Compressing $FILE"
	if command -v pigz >/dev/null 2>&1; then
	    pigz -p 4 -k "${FILE}"
	else
	    gzip -k "${FILE}"
	fi
    else
        info "Skipping $FILE (already compressed)"
    fi
done

# Run evaluation
header "Evaluating outputs..."
for FILE in ${OUTDIR}/*.{fasta,fastq}*; do
    B=$(basename "$FILE" | cut -f 1 -d.)
    N50=$(echo "$B" | cut -d _ -f 2)
    TOT=$(echo "$B" | cut -d _ -f 3)
    SIZE=$(echo "$B" | cut -d _ -f 4)

    EVAL=$(bin/n50 -b --csv "$FILE" | tail -n 1)
    TEST_N50=$(echo "$EVAL" | cut -f 4 -d,)
    TEST_TOT=$(echo "$EVAL" | cut -f 2 -d,)
    TEST_SIZE=$(echo "$EVAL" | cut -f 3 -d,)

    [[ "$N50" == "$TEST_N50" ]] && success "Correct N50 in $B" || fail "Wrong N50 in $B: expected $N50, got $TEST_N50"
    [[ "$TOT" == "$TEST_TOT" ]] && success "Correct total seqs in $B" || fail "Wrong total seqs in $B: expected $TOT, got $TEST_TOT"
    [[ "$SIZE" == "$TEST_SIZE" ]] && success "Correct total size in $B" || fail "Wrong total size in $B: expected $SIZE, got $TEST_SIZE"

    rm -f "$FILE"
done

rmdir "$OUTDIR" || true

