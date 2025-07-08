# `n50_simseqs` - Simulate Sequences

> A utility to generate simulated FASTA or FASTQ sequences with specified counts and lengths, useful for testing and benchmarking N50 calculation tools.

## Features

- Generates sequences in either FASTA or FASTQ format.
- Allows specification of multiple `COUNT*SIZE` pairs to create diverse datasets.
- Calculates and includes the N50 of the simulated sequences in the output filename.
- Supports custom output directory and filename prefix.

## Installation

### Prerequisites

- GCC compiler

### Compiling

To compile the program, use the following command:

```bash
gcc -o n50_simseqs src/n50_simseqs.c
```

## Usage

```bash
./n50_simseqs [--fasta|--fastq] -o OUTDIR [-p PREFIX] ARGS...
```

### Arguments

- `ARGS`: One or more specifications for the simulated reads in `COUNT*SIZE` format.
  - `COUNT`: The number of sequences to generate for this specification.
  - `SIZE`: The length of the sequences. Supports suffixes `K` (kilobases), `M` (megabases), `G` (gigabases).

### Options

- `--fasta`: Generate sequences in FASTA format (default).
- `--fastq`: Generate sequences in FASTQ format.
- `-o OUTDIR`: **Required**. The output directory where the generated file will be saved.
- `-p PREFIX`: Optional. A prefix to add to the output filename.
- `--version`: Show version number and exit.
- `-h`, `--help`: Show this help message and exit.

### Examples

1. Generate 1000 FASTA sequences of 500bp each into `output_dir`:

```bash
./n50_simseqs -o output_dir 1000*500
```

2. Generate a mix of FASTQ sequences: 5000 reads of 1Kbp and 1000 reads of 2.5Kbp, with a prefix `test_`:

```bash
./n50_simseqs --fastq -o my_sim_data -p test_ 5000*1K 1000*2500
```

3. Generate a single large FASTA file with 100 reads of 1Mbp:

```bash
./n50_simseqs -o large_data 100*1M
```

## Output

The program creates a single output file in the specified `OUTDIR`. The filename format is:

`[PREFIX][N50]_[TOTAL_READS]_[TOTAL_BASES].[fasta|fastq]`

- `N50`: The N50 value of the generated sequences.
- `TOTAL_READS`: The total number of sequences generated.
- `TOTAL_BASES`: The total number of bases generated.

Example filename: `test_500_6000_7500000.fastq`

## Version

`1.9.2`
