CC = gcc
CFLAGS = -Wall -Wextra -O3 
LDFLAGS = -lz -lpthread

SRC_DIR = src
BIN_DIR = bin
TEST_DIR = test
TARGET  = $(BIN_DIR)/n50 
COUNTBIN = $(BIN_DIR)/fqc
COUNTFABIN = $(BIN_DIR)/fac
COUNTFXBIN = $(BIN_DIR)/countfx
SIMTARGET = $(BIN_DIR)/gen
SIMDATA = test/sim/list.txt
# Find all n50 variant source files
N50_VARIANTS := $(wildcard $(SRC_DIR)/n50_*.c)
# Create target names for all n50 variants
N50_VARIANT_TARGETS := $(patsubst $(SRC_DIR)/n50_%.c,$(BIN_DIR)/n50_%,$(N50_VARIANTS))

.PHONY: all clean test

all: $(TARGET) $(SIMTARGET) $(TESTTARGET) $(N50_VARIANT_TARGETS) $(COUNTBIN) $(COUNTFABIN) $(COUNTFXBIN)

#Make targets
$(TARGET): $(SRC_DIR)/n50.c | $(BIN_DIR)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

$(TESTTARGET): $(SRC_DIR)/n50_opt.c | $(BIN_DIR)
	$(CC)  $(CFLAGS) $< -o $@ $(LDFLAGS)

$(SIMTARGET): $(SRC_DIR)/gen.c | $(BIN_DIR)
	$(CC) $(CFLAGS) $< -o $@  

$(COUNTBIN): $(SRC_DIR)/counts.c
	$(CC) -O3 -pthread $< -o $@ -lz -lpthread

$(COUNTFABIN): $(SRC_DIR)/countfa.c
	$(CC) -O3 -pthread $< -o $@ -lz -lpthread

$(COUNTFXBIN): $(SRC_DIR)/countfx.c
	$(CC) -O3 -pthread $< -o $@ -lz -lpthread

$(BIN_DIR):
	mkdir -p $(BIN_DIR)


# Rule for n50 variants
$(BIN_DIR)/n50_%: $(SRC_DIR)/n50_%.c | $(BIN_DIR)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

clean:
	rm -rf $(BIN_DIR)
	if [ -d "test/sim" ]; then echo "Removing sim"; rm -rf test/sim/; fi

$(SIMDATA): $(SIMTARGET)
	# Generate simulated files each with filename like {N50}_{num_seqs}_{total_length}.{format} \
	mkdir -p test/sim; \
	# ./program <min_seqs> <max_seqs> <min_len> <max_len> <tot_files> <format> <outdir>\
	# small datasets \
	$(SIMTARGET) 500        1000          5      2000000    5 fasta test/sim/; \
	$(SIMTARGET) 500        1000          5        10000    5 fastq test/sim/; \
	# large datasets
	$(SIMTARGET) 5000     10000          1000      20000000    3 fasta test/sim/; \
	$(SIMTARGET) 5000     10000          1000      20000000    3 fastq test/sim/; \
	for FILE in test/sim/*.{fasta,fastq}; do \
		gzip -k -f "$$FILE"; \
	done; \
	ls test/sim/*.* > $(SIMDATA)

# Test rule
test: $(TARGET) $(SIMTARGET)
	bash test/test.sh

# Original simple test
autotest: $(TARGET)
	@echo "Running simple test..."
	@echo ">seq1" > test.fasta
	@echo "ATCGATCGATCG" >> test.fasta
	@echo ">seq2" >> test.fasta
	@echo "ATCGATCGATCGATCGATCG" >> test.fasta
	@echo ">seq3" >> test.fasta
	@echo "ATCG" >> test.fasta
	@output=$$($(TARGET) test.fasta); \
	echo "Output: $$output"; \
	total_bp=$$(echo "$$output" | cut -f3); \
	total_seq=$$(echo "$$output" | cut -f4); \
	n50=$$(echo "$$output" | cut -f5); \
	if [ "$$total_bp" = "36" ] && [ "$$total_seq" = "3" ] && [ "$$n50" = "20" ]; then \
		echo "Simple test passed successfully!"; \
	else \
		echo "Simple test failed!"; \
		echo "Expected: 36 total bp, 3 sequences, N50 of 20"; \
		echo "Got: $$total_bp total bp, $$total_seq sequences, N50 of $$n50"; \
		exit 1; \
	fi
	@rm test.fasta
	@echo "Simple test completed."
	


benchmark: $(TARGET) $(SIMTARGET) $(SIMDATA)

	if [ -d "test/sim" ]; then \
		if [ ! -d "test/benchmark" ]; then mkdir -p test/benchmark; fi; \
		if [ ! -d "test/local/" ]; then mkdir -p test/local/; fi; \
		echo "Running simulation tests in test/sim"; \
		for file in test/sim/*.*; do \
			hyperfine --warmup 2 --max-runs 20 --export-csv "test/benchmark/$$(basename "$$file").csv" "$(TARGET) $$file" "seqfu stats $$file" "seqkit stats $$file"; \
		done; \
		hyperfine --warmup 2 --max-runs 20 --export-csv "test/benchmark/all.csv" "$(TARGET) test/sim/*.*" "seqfu stats test/sim/*.*" "seqkit stats test/sim/*.*"; \
		hyperfine --warmup 2 --max-runs 20 --export-csv "test/benchmark/all_fasta.csv" "$(TARGET) test/sim/*.fasta" "seqfu stats test/sim/*.fasta" "seqkit stats test/sim/*.fasta"; \
		hyperfine --warmup 2 --max-runs 20 --export-csv "test/benchmark/all_fastq.csv" "$(TARGET) test/sim/*.fastq" "seqfu stats test/sim/*.fastq" "seqkit stats test/sim/*.fastq"; \
		for FILE in test/sim/*.{fasta,fastq}; do \
			gzip -k -f "$$FILE"; \
		done; \
		hyperfine --warmup 2 --max-runs 20 --export-csv "test/benchmark/gz_all.csv" "$(TARGET) test/{sim,local}/*.gz" "seqfu stats test/{sim,local}/*.gz" "seqkit stats test/{sim,local}/*.gz"; \
		hyperfine --warmup 2 --max-runs 20 --export-csv "test/benchmark/gz_all_fasta.csv" "$(TARGET) test/{sim,local}/*.fasta.gz" "seqfu stats test/{sim,local}/*.fasta.gz" "seqkit stats test/{sim,local}/*.fasta.gz"; \
		hyperfine --warmup 2 --max-runs 20 --export-csv "test/benchmark/gz_all_fastq.csv" "$(TARGET) test/{sim,local}/*.fastq.gz" "seqfu stats test/{sim,local}/*.fastq.gz" "seqkit stats test/{sim,local}/*.fastq.gz"; \
	else \
		echo "test/sim directory does not exist. Skipping simulation tests."; \
	fi
 
