# `n50` - Calculate N50

> A fast and efficient tool for calculating N50 and other sequence statistics from FASTA and FASTQ files.

## Features


- Supports both FASTA and FASTQ formats
- Optimized for FASTQ raw file (Nanopore, PacBio)
- Handles gzipped input files

## Installation

### Prerequisites

- GCC compiler
- zlib library
- pthread library

### Compiling

To compile the program, use the following command:

```bash
make
```

or compile binaries like:

```bash
gcc -o n50 src/n50.c -lz -lpthread -O3
```

## Usage

```bash
./n50 [options] [filename]...
```

If no filename is provided, the program reads from standard input.

### Options

- `--fasta` or `-a`: Force FASTA input format
- `--fastq` or `-q`: Force FASTQ input format
- `--header` or `-H`: Print header in the output
- `--n50` or `-n`: Output only the N50 value
- `--version` and `--help` to display version and help information respectively

### Examples

1. Process a FASTA file:

```bash
./n50 input.fasta
```

2. Process a gzipped FASTQ file:

```bash
./n50 input.fastq.gz
```

3. Process a file with header output:
   
```bash
./n50 --header input.fasta
```

4. Output only the N50 value:

```bash
./n50 --n50 input.fasta
```

5. Process input from stdin:

```bash
cat input.fasta | ./n50 --format fasta
```

## Output

By default, the program outputs a tab-separated line with the following fields:

1. Format (FASTA or FASTQ)
2. Total sequence length
3. Total number of sequences
4. N50 value

When using the `--header` option, a header line is printed before the results.

When using the `--n50` option, only the N50 value is printed.

## Performance

The program uses multi-threading to process large files efficiently. It automatically adjusts the number of threads based on the input size, up to a maximum of 8 threads.

## Limitations

- The maximum number of threads is currently set to 8. This can be adjusted by modifying the `MAX_THREADS` constant in the source code.
- The initial capacity for storing sequence lengths is set to 1,000,000. For extremely large datasets, this value might need to be increased.
