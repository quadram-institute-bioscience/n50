import random

def gen_ctg_len(min_len, max_len, num_seqs):
    """
    Generate a list of contig lengths.
    """
    lengths = []
    for _ in range(num_seqs):
        lengths.append(random.randint(min_len, max_len))
    return lengths

def calculate_n50(lengths):
    """Calculate N50 from a list of sequence lengths."""
    sorted_lengths = sorted(lengths, reverse=True)
    total_length = sum(sorted_lengths)
    cumulative_length = 0
    for length in sorted_lengths:
        cumulative_length += length
        if cumulative_length >= total_length / 2:
            return (length, len(lengths), total_length)
        
def generate_contigs(N50, SUM_LEN, TOT_SEQS):
    if N50 > SUM_LEN or TOT_SEQS < 1:
        raise ValueError("Invalid input: N50 must be <= SUM_LEN and TOT_SEQS >= 1")
    
    # Step 2: Generate MAX_LEN, a random number between N50 and SUM_LEN
    MAX_LEN = random.randint(N50, SUM_LEN)
    print("MAX_LEN\t", MAX_LEN)
    
    contig_list = [MAX_LEN]  # First number in the list is MAX_LEN
    TMP_SUM = MAX_LEN        # Track the total sum of the contigs
    
    # Step 4: Add contigs to the list, decreasing in size
    while contig_list[-1] > N50:
        print(" TMP_SUM\t", TMP_SUM, "\t", TMP_SUM/SUM_LEN)
        # Generate a random contig smaller than the last one, but non-zero
        next_contig = random.randint(1, contig_list[-1])
        print("LAST_CTG\t", contig_list[-1])
        print(" NEW_CTG\t", next_contig, "\n")
        contig_list.append(next_contig)
        TMP_SUM += next_contig
        
    # Step 5: If TMP_SUM exceeds N50, replace the last added contig with N50
    if contig_list[-1] > N50:
        contig_list[-1] = N50
    
    # Step 6: Fill the remaining contigs until TMP_SUM == SUM_LEN
    while len(contig_list) < TOT_SEQS:
        remaining_sum = SUM_LEN - TMP_SUM
        if remaining_sum <= 0:
            break
        next_contig = min(remaining_sum, random.randint(1, MAX_LEN))
        contig_list.append(next_contig)
        TMP_SUM += next_contig

    return contig_list


def generate_sequence(length):
    """Generate a random DNA sequence of a given length."""
    return ''.join(random.choices('ACGT', k=length))


def write_fasta(sequences, outfile):
    """Write sequences in FASTA format."""
    try:
        with open(outfile, 'w') as f:
            for i, seq in enumerate(sequences):
                print(f">seq{i+1}", file=f)
                print(f"{seq}", file=f)
    except IOError as e:
        print("Unable to write {}:\n {}".format(outfile, e), file=sys.stderr)
    

def write_fastq(sequences, outfile):
    """Write sequences in FASTQ format (dummy quality scores)."""

    try:
        with open(outfile, 'w') as f:
            for i, seq in enumerate(sequences):
                print(f"@seq{i+1}")
                print(f"{seq}")
                print(f"+")
                print(f"{'I' * len(seq)}")
    except IOError as e:
        print("Unable to write {}:\n {}".format(outfile, e), file=sys.stderr)

def generate_sequences(lengths, format, outfile):
    """
    Generate sequences with a from a list of lengths.
    """
    # Generate sequences based on the generated lengths.
    sequences = [generate_sequence(length) for length in lengths]

    # Print the sequences to standard output in the chosen format.
    if format == 'FASTA':
        write_fasta(sequences, outfile)
    elif format == 'FASTQ':
        write_fastq(sequences, outfile)

if __name__ == "__main__":
    import argparse
    import sys, os
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-seqs", default=15, type=int, help="Minimum number of sequences")
    parser.add_argument("--max-seqs", default=5000, type=int, help="Maximum number of sequences")
    parser.add_argument("--min-len", default=12, type=int, help="Minimum length for a sequence")
    parser.add_argument("--max-len", default=100000,type=int, help="Maximum length for a sequence")
    parser.add_argument("--tot", default=10, type=int, help="Number of files to generate")
    parser.add_argument("--format", default="FASTA", help="Output format (FASTA or FASTQ)")
    parser.add_argument("--outdir", help="Output directory")
    args = parser.parse_args()

    if not args.outdir:
        print("Please provide an output directory", file=sys.stderr)
        # print help
        parser.print_help()
        sys.exit(1)
    if not os.path.exists(args.outdir):
        try:
            os.makedirs(args.outdir)
        except OSError as e:
            sys.stderr.write(str(e) + '\n')
            sys.exit(1)
    
    for i in range(args.tot):
        
        num_seqs = random.randint(args.min_seqs, args.max_seqs)
        contigs_lengths = gen_ctg_len(args.min_len, args.max_len, num_seqs)
        n50, num_seqs, total_length = calculate_n50(contigs_lengths)
        outfile = os.path.join(args.outdir, "{}_{}_{}.{}".format(n50, num_seqs, total_length, args.format.lower()))
        print("N50: {}\tNum_seqs: {}\tTot_len: {}".format(n50, num_seqs, total_length), file=sys.stderr)
        generate_sequences(contigs_lengths, args.format, outfile)
