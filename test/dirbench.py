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
            subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
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
                print(f"Adding {file} to {exts[ext]}", file=sys.stderr)
                
        
    return results

def filesize(filepath, unit="Mb"):
    """
    Return filesize in the required unit (b, kb, Mb, Gb).
    """
    size = os.path.getsize(filepath)
    if unit == "b":
        return size
    elif unit == "kb":
        return size / 1024
    elif unit == "Mb":
        return size / (1024 * 1024)
    elif unit == "Gb":
        return size / (1024 * 1024 * 1024)
    else:
        raise ValueError("Invalid unit. Choose from 'b', 'kb', 'Mb', 'Gb'.")

def gen_cmd(files, outdir, testname, params=['--max-runs', '5', '--warmup', '1']):
    cmd = ['hyperfine']
    csv = os.path.join(outdir, f"{testname}.csv")
    params.extend(['--export-csv', csv])
    cmd.extend(params)
    progs = ['n50', 'seqfu stats', 'seqkit stats --all']
    for p in progs:
        cmd.append(f"{p} {' '.join(files)}")
    
    return cmd

if __name__ == "__main__":
    import argparse
    args = argparse.ArgumentParser(description="Benchmark N50 calculation tools using hyperfine.")
    args.add_argument("FILES", nargs="+", help="N50 calculation tools to benchmark.")
    args.add_argument('-o','--outdir', help="")
    args.add_argument("--verbose", action="store_true", help="Print verbose output.")
    args = args.parse_args()
    if args.verbose:
        print(f"Received {len(args.FILES)} files", file=sys.stderr)
    check_deps(args.verbose)
    file_groups = split_files_type(args.FILES)
    for ftype in file_groups:
        print(ftype, "\t", len(file_groups[ftype]))
        if len(file_groups[ftype]) == 0:
            continue
        c = gen_cmd(file_groups[ftype], args.outdir, ftype)
        try:
            print(c)
            # run suppressing output
            subprocess.run(c, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except Exception as e:
            print(e)
            sys.exit(1)

    