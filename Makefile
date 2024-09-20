CC = gcc
CFLAGS = -Wall -Wextra -O3 
LDFLAGS = -lz -lpthread

SRC_DIR = src
BIN_DIR = bin
TEST_DIR = test
TARGET  = $(BIN_DIR)/n50 
SIMTARGET = $(BIN_DIR)/gen
SIMDATA = test/sim/list.txt
# Find all n50 variant source files
N50_VARIANTS := $(wildcard $(SRC_DIR)/n50_*.c)
# Create target names for all n50 variants
N50_VARIANT_TARGETS := $(patsubst $(SRC_DIR)/n50_%.c,$(BIN_DIR)/n50_%,$(N50_VARIANTS))

.PHONY: all clean test

all: $(TARGET) $(SIMTARGET) $(TESTTARGET) $(N50_VARIANT_TARGETS)

#Make targets
$(TARGET): $(SRC_DIR)/n50.c | $(BIN_DIR)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

$(TESTTARGET): $(SRC_DIR)/n50_opt.c | $(BIN_DIR)
	$(CC)  $(CFLAGS) $< -o $@ $(LDFLAGS)

$(SIMTARGET): $(SRC_DIR)/gen.c | $(BIN_DIR)
	$(CC) $(CFLAGS) $< -o $@  


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
	@passed=0; failed=0; \
	if [ -d "$(TEST_DIR)" ]; then \
		echo "[1] Running tests in $(TEST_DIR)"; \
		for file in $(TEST_DIR)/*.*; do \
			filename=$$(basename "$$file"); \
			expected_n50=$${filename%%.*}; \
			echo "Testing $$filename (Expected N50: $$expected_n50)"; \
			output=$$($(TARGET) "$$file"); \
			actual_n50=$$(echo "$$output" | cut -f 5); \
			if [ "$$actual_n50" = "$$expected_n50" ]; then \
				echo "  [OK] Got expected N50"; \
				passed=$$((passed + 1)); \
			else \
				echo "  FAIL: Expected N50: $$expected_n50, Got: $$actual_n50"; \
				failed=$$((failed + 1)); \
			fi; \
		done; \
	fi; \
	mkdir -p test/sim; \
	if [ -d "test/sim" ]; then \
		echo "[2] Generating simulated reads" \
		# Generate simulated files each with filename like {N50}_{num_seqs}_{total_length}.{format} \
		$(SIMTARGET) 10 35 1 2000 10 fasta test/sim/; \
		$(SIMTARGET) 5 15 1 1000 8 fastq test/sim/; \
		echo "Running simulation tests in test/sim"; \
		for file in test/sim/*.*; do \
			if echo "$$file" | grep -Eq '([0-9]+)_([0-9]+)_([0-9]+)\.(fasta|fastq)(\.gz)?$$'; then \
				n50=$$(echo "$$file" | sed -E 's#.*/([0-9]+)_([0-9]+)_([0-9]+)\..*#\1#'); \
				seqs=$$(echo "$$file" | sed -E 's#.*/([0-9]+)_([0-9]+)_([0-9]+)\..*#\2#'); \
				totlen=$$(echo "$$file" | sed -E 's#.*/([0-9]+)_([0-9]+)_([0-9]+)\..*#\3#'); \
				ext=$$(echo "$$file" | sed -E 's/.*\.([^.]+)(\.gz)?$$/\1/'); \
				gz=$$(echo "$$file" | sed -E 's/.*(\.(gz))?$$/\1/'); \
				echo "Testing $$(basename "$$file") (N50: $$n50, Seqs: $$seqs, TotLen: $$totlen)"; \
				output=$$($(TARGET) --"$$ext" "$$file"); \
				actual_size=$$(echo "$$output" | cut -f 3); \
				actual_n50=$$(echo "$$output" | cut -f 5); \
				actual_seqs=$$(echo "$$output" | cut -f 4); \
				if [ "$$actual_seqs" -ne "$$seqs" ]; then \
					echo "  [FAIL]    Expected $$seqs sequences, but found $$actual_seqs"; \
					failed=$$((failed + 1)); \
				else \
					echo "  [OK]      Sequence count is $$seqs"; \
					passed=$$((passed + 1)); \
				fi; \
				if [ "$$actual_n50" -ne "$$n50" ]; then \
					echo "  [FAIL]    Expected N50 of $$n50, but found $$actual_n50"; \
					failed=$$((failed + 1)); \
				else \
					echo "  [OK]      N50 is $$n50"; \
					passed=$$((passed + 1)); \
				fi; \
				if [ "$$actual_size" -ne "$$totlen" ]; then \
					echo "  [FAIL]    Expected size of $$actual_size, but found $$totlen"; \
					failed=$$((failed + 1)); \
				else \
					echo "  [OK]      Total bp $$actual_size"; \
					passed=$$((passed + 1)); \
				fi; \
			fi; \
		done; \
		rm test/sim/*; \
	else \
		echo "test/sim directory does not exist. Skipping simulation tests."; \
	fi; \
	echo "Tested $(TARGET)"; \
	echo "Tests completed. Passed: $$passed, Failed: $$failed"; \
	if [ $$failed -ne 0 ]; then exit 1; fi

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
 