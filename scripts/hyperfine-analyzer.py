import pandas as pd
import glob
from pathlib import Path
"""
    Analyzes benchmark results from multiple hyperfine CSV output files.
    
    This function processes multiple CSV files containing hyperfine benchmark results,
    calculating aggregate statistics and determining which command was fastest across
    all runs. It provides two main analyses:
    1. The mean execution time for each command across all benchmark files
    2. A count of how many times each command was the fastest in individual benchmarks
    
    Args:
        file_pattern (str): Glob pattern to match hyperfine CSV output files 
                           (e.g., "*.csv" or "benchmarks/*.csv")
    
    Returns:
        None - Results are printed to stdout
        
    Raises:
        FileNotFoundError: If no files match the provided pattern
        
    Example:
        >>> analyze_hyperfine_results("benchmarks/*.csv")
        Aggregate mean execution times:
        command1: 1.234 seconds
        command2: 2.345 seconds
        
        Number of times each command was fastest:
        command1: 3 times
        command2: 1 times
"""
def analyze_hyperfine_results(file_pattern):
    # Read all CSV files matching the pattern
    all_data = []
    for file in glob.glob(file_pattern):
        df = pd.read_csv(file)
        all_data.append(df)
    
    # Combine all dataframes
    combined_df = pd.concat(all_data, ignore_index=True)
    
    # 1. Calculate aggregate mean execution time for each command
    aggregate_means = combined_df.groupby('command')['mean'].mean().round(3)
    
    # 2. Count how many times each command was the fastest
    wins_count = {}
    for file in glob.glob(file_pattern):
        df = pd.read_csv(file)
        fastest = df.loc[df['mean'].idxmin(), 'command']
        wins_count[fastest] = wins_count.get(fastest, 0) + 1
    
    # Print results
    print("Aggregate mean execution times:")
    for command, mean in aggregate_means.items():
        print(f"{command}: {mean:.3f} seconds")
    
    print("\nNumber of times each command was fastest:")
    for command in aggregate_means.index:
        wins = wins_count.get(command, 0)
        print(f"{command}: {wins} times")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='Analyze multiple hyperfine benchmark CSV files'
    )
    parser.add_argument(
        'pattern', 
        help='File pattern for CSV files (e.g., "*.csv" or "benchmarks/*.csv")'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Print additional details'
    )

    args = parser.parse_args()

    try:
        # Check if any files match the pattern
        matching_files = list(glob.glob(args.pattern))
        if not matching_files:
            raise FileNotFoundError(f"No files found matching pattern: {args.pattern}")
        
        if args.verbose:
            print(f"Processing {len(matching_files)} files:")
            for f in matching_files:
                print(f"  - {f}")
            print()

        analyze_hyperfine_results(args.pattern)

    except Exception as e:
        print(f"Error: {str(e)}")
        exit(1)

