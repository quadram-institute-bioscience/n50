#!/usr/bin/env python
"""
A program to benchmark N50 calculations tools 
using hyperfine.
"""

import os
import sys
import subprocess

def check_deps(verbose):
    """
    Check if hyperfine is installed.
    """
    test_cmds = [
        ['hyperfine', '--version'],
        ['seqfu', '--version'],
        ['n50', '--version'],
        ['seqkit', 'version'],
    ]

    for cmd in test_cmds:
        try:
            subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            if verbose:
                print(f"Command {cmd[0]} found.")
        except FileNotFoundError:
            print(f"Command {cmd[0]} not found. Exiting.")
            sys.exit(1)

def split_files_type(files):
    """
    Split the list of files into a dictionary of lists
    """

    results = {
        'fasta' : [],
        'fastq' : [],
        'fasta_gz' : [],
        'fastq_gz' : [],
    }
    exts = {
        '.fasta' : 'fasta',
        '.fastq' : 'fastq',
        '.fasta.gz' : 'fasta_gz',
        '.fastq.gz' : 'fastq_gz',
        '.fq' : 'fastq',
        '.fq.gz' : 'fastq_gz',
        '.fa' : 'fasta',
        '.fa.gz' : 'fasta_gz',
        '.fna'   : 'fasta'
    }
    for file in files:
        for ext in exts:
            if file.endswith(ext):
                results[exts[ext]].append(file)
        
    return results

def gen_cmd(files, params=['--max-runs 5', '--warmup 1']):
    cmd = ['hyperfine']
    cmd.extend(params)
    progs = ['n50', 'seqfu stats', 'seqkit stats --all']
    for p in progs:
        cmd.append(f"{p} {' '.join(files)}")
    
    print(cmd)

if __name__ == "__main__":
    import argparse
    args = argparse.ArgumentParser(description="Benchmark N50 calculation tools using hyperfine.")
    args.add_argument("FILES", nargs="+", help="N50 calculation tools to benchmark.")
    args.add_argument('-o','--outdir', help="")
    args.add_argument("--verbose", action="store_true", help="Print verbose output.")
    args = args.parse_args()
    check_deps(args.verbose)
    file_groups = split_files_type(args.FILES)
    for ftype in file_groups:
        print(ftype, "\t", len(file_groups[ftype]))
        c = gen_cmd(file_groups[ftype])

    