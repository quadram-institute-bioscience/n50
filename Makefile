CC = gcc
CFLAGS = -Wall -Wextra -O3
LDFLAGS = -lz -lpthread

SRC_DIR = src
BIN_DIR = bin
TEST_DIR = test
TARGET = $(BIN_DIR)/n50

.PHONY: all clean test

all: $(TARGET)

$(TARGET): $(SRC_DIR)/n50.c | $(BIN_DIR)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

clean:
	rm -rf $(BIN_DIR)

test: $(TARGET)
	@echo "Running tests..."
	@passed=0; failed=0; \
	for file in $(TEST_DIR)/*.*; do \
		filename=$$(basename $$file); \
		expected_n50=$${filename%%.*}; \
		echo "Testing $$filename (Expected N50: $$expected_n50)"; \
		output=$$($(TARGET) $$file); \
		actual_n50=$$(echo "$$output" | cut -f4); \
		if [ "$$actual_n50" = "$$expected_n50" ]; then \
			echo "  Passed"; \
			passed=$$((passed + 1)); \
		else \
			echo "  Failed. Expected N50: $$expected_n50, Got: $$actual_n50"; \
			failed=$$((failed + 1)); \
		fi; \
	done; \
	echo "Tests completed. Passed: $$passed, Failed: $$failed"; \
	if [ $$failed -ne 0 ]; then exit 1; fi

# Original simple test
test-simple: $(TARGET)
	@echo "Running simple test..."
	@echo ">seq1" > test.fasta
	@echo "ATCGATCGATCG" >> test.fasta
	@echo ">seq2" >> test.fasta
	@echo "ATCGATCGATCGATCGATCG" >> test.fasta
	@echo ">seq3" >> test.fasta
	@echo "ATCG" >> test.fasta
	@output=$$($(TARGET) test.fasta); \
	echo "Output: $$output"; \
	total_bp=$$(echo "$$output" | cut -f2); \
	total_seq=$$(echo "$$output" | cut -f3); \
	n50=$$(echo "$$output" | cut -f4); \
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
