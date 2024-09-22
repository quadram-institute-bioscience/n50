# `gen`: DNA Sequence Generator

## Overview

The DNA Sequence Generator is a C program designed to create random DNA sequences in FASTA or FASTQ format focusing on some **parameters focusing on their lengths**.
It generates multiple files containing sequences of varying lengths.

This tool is particularly useful for:

- Testing tools to gather statistics from DNA sequences
- Benchmarking bioinformatics pipelines

## Features

- Calculates and includes N50 statistics in the output filenames
- Supports both FASTA and FASTQ output formats
- Deterministic output based on a provided seed for reproducibility
- Output filename format: `N50_TOTSEQS_SUMLEN.{fasta|fastq}` to make easy to test the N50 calculation

## Compilation

To compile the program, use a C compiler such as GCC:

```bash
gcc -o dna_generator dna_generator.c -lm
```

## Usage

```bash
./gen <min_seqs> <max_seqs> <min_len> <max_len> <tot_files> <format> <outdir> 
```

### Parameters

- `<min_seqs>`: Minimum number of sequences per file
- `<max_seqs>`: Maximum number of sequences per file
- `<min_len>`: Minimum length of each sequence
- `<max_len>`: Maximum length of each sequence
- `<tot_files>`: Total number of files to generate
- `<format>`: Output format (either "fasta" or "fastq")
- `<outdir>`: Directory to store the output files
- `<seed>`: Seed for the random number generator (for reproducibility)

### Example

```bash
gen  10 100 1000 10000 5 fasta output_dir 
```

This command will generate 5 FASTA files in the `output_dir` directory. Each file will contain between 10 and 100 sequences, with lengths ranging from 1000 to 10000 base pairs. The random number generator will be initialized with a static seed to ensure reproducibility.

## Output

The program generates files with names in the format:

```text
${N50}_${TOTSEQS}_${SUMLEN}.{fasta|fastq}
```

Where:

- `N50` is the N50 statistic of the sequences in the file
- `TOTSEQS` is the total number of sequences in the file
- `SUMLEN` is the sum of all sequence lengths in the file

### FASTA Format

In FASTA format, each sequence is represented as:

```text
>seq1
ATCGATCGATCG...
>seq2
GCTAGCTAGCTA...
```

### FASTQ Format

In FASTQ format, each sequence is represented as:

```text
@seq1
ATCGATCGATCG...
+
IIIIIIIIIIII...
@seq2
GCTAGCTAGCTA...
+
IIIIIIIIIIII...
```

Note: The quality scores in FASTQ format are dummy values (all 'I') for simplicity.

## License

This program is provided under the MIT License. See the source code for full license text.

## Author

Andrea Telatin, 2023
Quadram Institute Bioscience

