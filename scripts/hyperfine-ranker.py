import pandas as pd
import sys
from collections import defaultdict

def analyze_hyperfine_results(files):
    # Read all provided CSV files
    all_data = []
    for file in files:
        try:
            df = pd.read_csv(file)
            all_data.append(df)
        except Exception as e:
            print(f"Error reading {file}: {e}", file=sys.stderr)
            continue
    
    if not all_data:
        print("No valid data to process.")
        return
    
    # Combine all dataframes
    combined_df = pd.concat(all_data, ignore_index=True)
    
    # 1. Calculate aggregate mean execution time for each command
    aggregate_means = combined_df.groupby('command')['mean'].mean().round(3)
    
    # 2. Count rankings for each command
    rankings = defaultdict(lambda: defaultdict(int))
    for df in all_data:
        if 'mean' not in df.columns or 'command' not in df.columns:
            continue
        # Sort by mean execution time and get rankings
        ranked = df.sort_values('mean')['command'].tolist()
        for position, command in enumerate(ranked, 1):
            rankings[command][position] += 1
    
    # Print results
    print("Aggregate mean execution times:")
    for command, mean in aggregate_means.items():
        print(f"{command}: {mean:.3f} seconds")
    
    print("\nRankings distribution:")
    # Get total number of commands to know the possible ranks
    max_rank = len(aggregate_means)
    
    for command in aggregate_means.index:
        rank_counts = [rankings[command][i] for i in range(1, max_rank + 1)]
        rank_str = " | ".join(f"#{i}: {count}" for i, count in enumerate(rank_counts, 1))
        print(f"{command}: {rank_str}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='Analyze multiple hyperfine benchmark CSV files'
    )
    parser.add_argument(
        'files', 
        nargs='+',
        help='Paths to CSV files (e.g., benchmarks/file1.csv benchmarks/file2.csv)'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Print additional details'
    )

    args = parser.parse_args()

    try:
        if args.verbose:
            print(f"Processing {len(args.files)} files:")
            for f in args.files:
                print(f"  - {f}")
            print()

        analyze_hyperfine_results(args.files)

    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)

