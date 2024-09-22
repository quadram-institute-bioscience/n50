# n50_binner

This program analyzes a FASTQ file and counts the number of reads falling into predefined length bins. 

It can be used to generate a simplified summary of the read length distribution in a 
FASTQ file, to be then used to 
[generate a similar FASTX file](README_N50_GENERATE.md).

## What it does

The program:

1. Reads a (compressed) FASTQ file.
2. Counts the number of reads falling into each of 16 predefined length bins.
3. Outputs a CSV-formatted result showing the number of reads in each bin.

## How to compile

To compile the program, you need a C compiler (such as gcc) and the zlib library installed. Use the following command:

```bash
gcc -o fastq_length_counter fastq_length_counter.c -lz
```

This will create an executable named `fastq_length_counter`.

## How to use

Run the program from the command line, providing the path to a gzip-compressed FASTQ file as an argument:

```bash
./n50_binner path/to/your/file.fastq.gz
```

The program will process the file and output the results to stdout in CSV format.

## Examples

1. Running the program:

```bash
./fastq_length_counter sample.fastq.gz
```

2. Saving the output to a file:

```bash
./fastq_length_counter sample.fastq.gz > length_distribution.csv
```

3. Processing multiple files:

```bash
for file in *.fastq.gz; do
    echo "Processing $file"
    ./fastq_length_counter "$file" > "${file%.fastq.gz}_length_distribution.csv"
done
```

## Output format

The output is in CSV format with two columns:
1. Bin: The upper limit of the length bin
2. Number of Reads: The count of reads falling into that bin

Example output:

```text
Bin,Number of Reads
10,0
100,5
1000,1000
2500,5000
...
```

## Notes

- The program uses 16 predefined bins: 10, 100, 1000, 2500, 5000, 10000, 20000, 35000, 50000, 75000, 100000, 200000, 300000, 500000, 750000, and 1000000.
- Reads longer than 1,000,000 bases will be counted in the last bin.
- The program assumes the input file is in the standard FASTQ format and is gzip-compressed.

## Dependencies

- zlib library for reading gzip-compressed files

Remember to install zlib development files before compiling. On Ubuntu or Debian, you can do this with:

```bash
sudo apt-get install zlib1g-dev
```

## License

This program is provided under the MIT License. See the source code for full license text.

## Author

Andrea Telatin, 2023
Quadram Institute Bioscience
