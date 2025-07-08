# `n50` - Calculate N50

> A fast and efficient tool for calculating N50 and other sequence statistics from FASTA and FASTQ files.

## Features

- Supports both FASTA and FASTQ formats (auto-detected)
- Handles gzipped input files
- Multi-threaded processing for efficiency
- Outputs in TSV, CSV, or JSON format
- Calculates N50, N75, N90, I50, GC content, total sequences, total length, average, minimum, and maximum sequence lengths.

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
./n50 [options] FILES...
```

If no filename is provided, the program reads from standard input (`-`).

### Arguments

- `FILES`: One or more FASTA/FASTQ files, gzipped or not. Use `-` for STDIN.

### Options

- `-a`, `--abs`: Print file paths as absolute paths in the output.
- `-b`, `--basename`: Print file paths as basename only (e.g., `file.fq.gz`) in the output.
- `-j`, `--json`: Output results in JSON format.
- `-c`, `--csv`: Output results in CSV format (default is TSV).
- `-h`, `--help`: Show this help message and exit.
- `-v`, `--version`: Show version number and exit.

### Examples

1. Process a FASTA file:

```bash
./n50 input.fasta
```

2. Process a gzipped FASTQ file and output in CSV format:

```bash
./n50 -c input.fastq.gz
```

3. Process multiple files and output in JSON format:
   
```bash
./n50 -j file1.fasta file2.fastq.gz
```

4. Process input from stdin:

```bash
cat input.fasta | ./n50 -
```

## Output

By default, the program outputs a tab-separated line (TSV) for each input file.

### TSV/CSV Output Columns:

1.  `Filepath`: Path to the input file.
2.  `TotSeqs`: Total number of sequences.
3.  `TotLen`: Total length of all sequences.
4.  `N50`: N50 value.
5.  `N75`: N75 value.
6.  `N90`: N90 value.
7.  `I50`: Number of sequences that contribute to N50.
8.  `GC`: GC content percentage.
9.  `Avg`: Average sequence length.
10. `Min`: Minimum sequence length.
11. `Max`: Maximum sequence length.

### JSON Output Format:

When using the `--json` option, the output is an array of JSON objects, where each object represents the statistics for a file. The keys are: `File`, `TotSeqs`, `TotLen`, `N50`, `N75`, `N90`, `I50`, `GC`, `Avg`, `Min`, `Max`.

## Performance

The program uses multi-threading to process large files efficiently. It utilizes up to 4 threads by default.

## Version

`1.9.2`