# sim-reads

sim-reads is a command-line tool for generating simulated genomic read datasets. 
It supports both Illumina-like and Nanopore-like data generation, allowing researchers to create datasets tailored to their analysis needs.

Will do:

- **Illumina-like Datasets:** Generate reads with fixed lengths (100, 150, 200, 300 bp).
- **Nanopore-like Datasets:** Generate reads with variable lengths based on predefined proportions.
- **Compression:** Automatically compresses output FASTA and FASTQ files using gzip.
- **Customizable:** Specify output directories, binary paths, and total dataset size.

## Requirements

- Python 3.6+
- `n50_simreads2` binary
- Required Python packages:
  - `subprocess`
  - `os`
  - `typing`
  - `gzip`
  - `shutil`
  - `argparse`
 
## Usage

Generate simulated reads with default settings:bash
python generate_simulated_reads.py

Specify a custom output directory and dataset size:bash
python generate_simulated_reads.py -o data/sim_reads -m 500

Specify a custom path to the `n50_simreads2` binary

```bash
python sim-reads.py -b /usr/local/bin/n50_simreads2
```
