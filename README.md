# krakenXtract

[![Release](https://github.com/Sam-Sims/krakenXtract/actions/workflows/release.yaml/badge.svg)](https://github.com/Sam-Sims/krakenXtract/actions/workflows/release.yaml)
![GitHub release (with filter)](https://img.shields.io/github/v/release/sam-sims/krakenxtract)
![crates.io](https://img.shields.io/crates/v/krakenxtract
)

Extract reads from a FASTQ file based on taxonomic classification via Kraken2.

Written in Rust.

## Motivation

Heavily inspired by the great [KrakenTools](https://github.com/jenniferlu717/KrakenTools). 

Having been wanting to experiment with Rust for a while, this is essentially an implementation of the `extract_kraken_reads.py` script, [re-written](https://www.reddit.com/media?url=https%3A%2F%2Fi.redd.it%2Fgood-for-you-crab-v0-5v9ygeh9r1c91.jpg%3Fs%3Dd759db5275e32c6e2bd5c22bddbd783acca46247) in Rust. The main motivation was to provide a speedup when extracting a large number of reads from large FASTQ files - and to learn Rust!

This is currently an early implementation, with plans to expand functionality.

## Current features

- Extract all reads from a `fastq` file based on a taxonomic id
- Extract all the parents or the children of the specified taxon id
- Supports interleaved `fastq` files
- Supports both uncompressed or `gzip` inputs.
- Multithreaded
- ~4.4x speed up compared to KrakenTools

### Benchmarks (WIP)

For more detail see [benchmarks](benchmarks/benchmarks.md)

## Installation

Download the latest release.

Alternatively, build from source:

**Clone the repository:**

```bash
git clone https://github.com/Sam-Sims/krakenxtract
```

**Install rust/cargo:**

To install please refer to the rust documentation: [docs](https://doc.rust-lang.org/cargo/getting-started/installation.html)

**Build and add to path:**

```bash
cd kraken-extract
cargo build --release
export PATH=$PATH:$(pwd)/target/release
```

All executables will be in the directory kraken-extract/target/release.

## Usage

```bash
krakenXtract --kraken <kraken_output> --fastq <fastq_file> --taxid <taxonomic_id> --output <output_file>
```

## Arguments

```
-k, --kraken <KRAKEN_OUTPUT>            
-t, --taxid <TAXID>              
-r, --report <REPORT_OUTPUT>            
-f, --fastq <FASTQ_FILE>              
-o, --output <OUTPUT_LOCATION>            
--compression-mode <COMPRESSION>      [default: fast]         
--parents                    
--children                   
--no-compress                
--exclude                    
  -h, --help                       Print help
  -V, --version                    Print version
```

`--parents`: This will extract all the reads classified at all taxons between the root and the specified `--taxid`

`--children`: This will extract all the reads classified as decendents or subtaxa of `--taxid` (Including the taxid)

`--compression_mode`: This defines the compression mode of the output `fastq.gz` file - fast / default / best

`--no-compress`: This will output a plaintext `fastq` file

`--exclude`: This will output every read except those matching the taxid. Works with `--parents` and `--children`

## Future plans

- [x] Support unzipped fastq files
- [ ] Support paired end FASTQ files
- [x] `--include-parents` and `--include-children` arguments
- [ ] Supply multiple taxonomic IDs to extract
- [x] Exclude taxonomic IDs
- [ ] `--append`
- [x] `--compression-mode <fast/default/best>`
- [x] More verbose output
- [ ] Proper benchmarks
- [ ] Output fasta format (for blast??)
- [x] Output non `gz`
- [ ] Tests

## Version

- 0.2.3

## Changelog

### 0.2.3

- Code optimisations

### 0.2.2

- Increased verbosity of outputs to user
- `--no-compress` flag to output a standard, plaintext fastq file
- `--exclude` to exclude specified reads. Works with `--children` and `--parents`
- Docstrings under the hood

### 0.2.1

- Fixes to reduce memory usage

### 0.2.0

- Detect and handle `gz` files or plain files
- `--compression` arg to select compression type
- `zlib-ng` to speed up gzip handling
- `--children` and `--parents` to save children and parents based on kraken report

### 0.1.0

- First release
