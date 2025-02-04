# hyperfine-analyzer

A command-line tool for analyzing multiple [hyperfine](https://github.com/sharkdp/hyperfine) benchmark results, providing aggregate statistics and determining which commands perform best across multiple benchmark runs.
 
- Calculates aggregate mean execution times across multiple benchmark files
- Tracks how often each command was the fastest
- Supports glob patterns for processing multiple benchmark files
- Simple command-line interface with verbose output option

## Usage

```bash
python hyperfine-analyzer.py "benchmarks/*.csv"
python hyperfine-analyzer.py -v "*.csv"
```

### Requirements
* Python 3.6+
* pandas
