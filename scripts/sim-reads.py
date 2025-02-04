#!/usr/bin/env python3
import subprocess
import os
from typing import List, Tuple
import gzip
import shutil
import argparse

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Generate simulated read datasets')
    parser.add_argument('-o', '--outdir', default='simulated_reads',
                      help='Output directory (default: simulated_reads)')
    parser.add_argument('-b', '--binary', default='bin/n50_simreads2',
                      help='Path to n50_simreads binary (default: bin/n50_simreads)')
    parser.add_argument('-m', '--megabases', type=int, default=200,
                      help='Approximate total size in megabases (default: 200)')
    return parser.parse_args()

def calculate_sequences(target_size: int, read_length: int) -> int:
    """Calculate number of sequences needed to reach target size."""
    return target_size // read_length

def generate_illumina_params(target_size: int) -> List[str]:
    """Generate parameters for Illumina-like datasets."""
    read_lengths = [100, 150, 200, 300]
    return [f"{calculate_sequences(target_size, length)}*{length}" 
            for length in read_lengths]

def generate_nanopore_params(target_size: int) -> List[str]:
    """Generate parameters for Nanopore-like datasets with variable read lengths."""
    length_ranges = [

        (      1_000,     10_000, 0.4),  # 40% short reads
        (     10_001,    100_000, 0.3),  # 40% medium reads
        (    100_001,  1_000_000, 0.1),  # 10% long reads
        (  1_000_001,  2_000_000, 0.1),  # 10% long reads
        (  2_000_001, 50_000_000, 0.1)   # 10% long reads
    ]
    
    params = []
    remaining_size = target_size
    
    for min_len, max_len, proportion in length_ranges:
        avg_len = (min_len + max_len) // 2
        range_size = int(target_size * proportion)
        num_seqs = range_size // avg_len
        params.append(f"{num_seqs}*{avg_len}")
        remaining_size -= num_seqs * avg_len
    
    return params

def run_tool(binary: str, params: List[str], output_dir: str, prefix: str, 
             file_format: str) -> Tuple[str, str]:
    """Run the simreads tool with given parameters."""
    cmd = [binary, f"--{file_format}", "-o", output_dir, 
           "-p", prefix] + params
    print(f"Running command: {' '.join(cmd)}")
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    
    # Parse the output to find the actual output file
    output_line = [line for line in result.stdout.split('\n') 
                  if line.startswith('Output written to:')]
    if not output_line:
        raise RuntimeError("Could not find output file in tool output")
    
    output_file = output_line[0].split(': ')[1].strip()
    gz_output_file = f"{output_file}.gz"
    return output_file, gz_output_file

def gzip_file(input_file: str, output_file: str):
    """Gzip the input file."""
    with open(input_file, 'rb') as f_in:
        with gzip.open(output_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

def main():
    args = parse_args()
    
    # Convert megabases to bases
    target_size = args.megabases * 1_000_000
    
    # Create output directory
    os.makedirs(args.outdir, exist_ok=True)
    
    # Generate Illumina-like datasets
    for read_length in [100, 150, 200, 300]:
        prefix = f"illumina_{read_length}_"  # Added trailing underscore
        params = [f"{calculate_sequences(target_size, read_length)}*{read_length}"]
        
        # Generate FASTA
        fasta_file, fasta_gz = run_tool(args.binary, params, args.outdir, prefix, "fasta")
        gzip_file(fasta_file, fasta_gz)
        
        # Generate FASTQ
        fastq_file, fastq_gz = run_tool(args.binary, params, args.outdir, prefix, "fastq")
        gzip_file(fastq_file, fastq_gz)
    
    # Generate Nanopore-like datasets
    nanopore_params = generate_nanopore_params(target_size)
    prefix = "nanopore_mixed_"  # Added trailing underscore
    
    # Generate FASTA
    fasta_file, fasta_gz = run_tool(args.binary, nanopore_params, args.outdir, prefix, "fasta")
    gzip_file(fasta_file, fasta_gz)
    
    # Generate FASTQ
    fastq_file, fastq_gz = run_tool(args.binary, nanopore_params, args.outdir, prefix, "fastq")
    gzip_file(fastq_file, fastq_gz)

if __name__ == "__main__":
    main()
