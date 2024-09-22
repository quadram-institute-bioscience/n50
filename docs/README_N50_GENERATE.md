# n50_generate

This program processes an input file containing read length distribution data and prepares it for use with [n50_simreads](README_N50_SIMREADS.md). 
It generates a string representation of the read length distribution and executes n50_simreads with the prepared data.

## Features

- Processes input files with read length distribution data from [n50_binner](README_N50_BINNER.md)
- Runs [n50_simreads](README_N50_SIMREADS.md) to generate reads based on the input data

## Usage

```bash
n50_prepare -i INPUTFILE -o OUTDIR [-f FORMAT] [-s PATH] [-v]
```

### Options:

- `-i INPUTFILE` : Path to the input file (required)
- `-o OUTDIR`    : Output directory for n50_simreads results (required)
- `-f FORMAT`    : Output format (optional, FASTQ by default, FASTA also supported)
- `-s PATH`      : Path to n50_simreads executable (optional)
- `-v`           : Verbose mode, prints additional information
- `-h`           : Display help message

## Input File Format

The input file should be a CSV file with the following format:

```text
length,count
100,1000
200,500
300,250
...
```

The first line (header) is skipped during processing.

## How to Compile

To compile the program, use a C compiler such as gcc:

```bash
gcc -o n50_prepare n50_prepare.c
```

This will create an executable named `n50_prepare`.

## Examples

1. Basic usage:

```bash
./n50_prepare -i input_distribution.csv -o output_directory
```

1. Specifying FASTA output format:

```bash
./n50_prepare -i input_distribution.csv -o output_directory -f FASTA
```

1. Using a custom path for n50_simreads:

```bash
./n50_prepare -i input_distribution.csv -o output_directory -s /path/to/n50_simreads
```

1. Running in verbose mode:

```bash
./n50_prepare -i input_distribution.csv -o output_directory -v
```

## Output

The n50_simreads output will be saved in the specified output directory.

## Statistics

The program calculates and displays:

- Total number of reads
- Maximum read length

## Dependencies

- n50_simreads (should be in the same directory as n50_prepare or specified with -s option)

## Notes

- The program assumes that n50_simreads is in the same directory as n50_prepare unless specified otherwise, but you can supply a custom path with the `-s` option.
- Make sure you have the necessary permissions to execute `n50_simreads` and write to the output directory.
- Invalid data in the input file will be skipped with a warning message.

## License

This program is provided under the MIT License. See the source code for full license text.

## Author

Andrea Telatin, 2023
Quadram Institute Bioscience
