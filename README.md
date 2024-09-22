# N50 Calculator v2.0

[![N50 Calculator CI](https://github.com/quadram-institute-bioscience/n50/actions/workflows/test_n50.yml/badge.svg)](https://github.com/quadram-institute-bioscience/n50/actions/workflows/test_n50.yml)

A fast and efficient tool for calculating N50 and other sequence statistics from FASTA and FASTQ files.

## Tools

- [`n50` calculate N50](docs/README_N50.md)

- [`n50_simreads`](docs/README_N50_SIMREADS.md), to simulate reads based on the lengths desired
- [`n50_binner`](docs/README_N50_BINNER.md),  to generate a summary of reads lengths from a FASTQ file to be used with `n50_generate`
- [`n50_generate`](docs/README_N50_GENERATE.md) uses the output of `n50_binner` to generate reads (using `n50_simreads`) 
  
- [`gen`](docs/README_GEN.md), alternative generator
- [benchmark notes](docs/README_BENCHMARK.md)


## General requirements

- C compiler
- zlib, pthread libraries

## Compiling

```bash
make all
make test
```

## Author

Andrea Telatin, 2023

## License

This program is open-source software released under the MIT License

## Contributing

Contributions to improve the N50 Calculator are welcome.
Please submit pull requests or open issues on the project's repository.
Be kind and adhere to the Code of Conduct  like the [Contributor Covenant](https://www.contributor-covenant.org/).
