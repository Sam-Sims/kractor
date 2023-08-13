# kraken-extract

Extract reads from a FASTQ file based on taxonomic classification via Kraken2.

Written in Rust.

## Background

I recently wanted to extract reads from a large-ish (6GB) FASTQ file (~5.5 million reads), based on taxonomic classifications. For that I used the great [KrakenTools](https://github.com/jenniferlu717/KrakenTools). This however took a while both parse the Kraken2 output file and extract/write the matching reads. Having been wanting to experiment with Rust for a while, this inspired me to re-implement the `extract_kraken_reads.py` script in Rust as a learning exercise.

This is currently an early implementation (and my first Rust programme!), with plans to expand functionality.

## Current features

- Extract all reads from a `fastq` file based on a taxonomic id
- Extract parents or the children of the specified taxon id
- Supports both uncompressed or `gzip` inputs.
- Multithreaded
- ~ 300% speed up over KrakenTools  

### Benchmarks (rough)

Based on 6.1Gb fastq.gz with 5,454,495 reads | 1.8Gb kraken output - 

Time to parse the kraken output, extract all matching reads, and write to new fastq file.

**KrakenTools (Output non gzip):**
| Type | Time       |
|------|------------|
| real | 6m 49s |

**kraken-extract:**
| Type | Time    |
|------|---------|
| real | 2m 04s |

## Installation

Clone the repository:

```bash
git clone https://github.com/Sam-Sims/kraken-extract
```

Build from source:

```bash
cd kraken-extract
cargo build --release
export PATH=$PATH:$(pwd)/target/release
```

All executables will be in the directory kraken-extract/target/release.

## Usage

```bash
kraken-extract --kraken <kraken_output> --fastq <fastq_file> --taxid <taxonomic_id> --output <output_file>
```

## Arguments

```
-k, --kraken <KRAKEN_OUTPUT>            
-t, --taxid <TAXID>              
-r, --report <REPORT_OUTPUT>            
-f, --fastq <FASTQ_FILE>              
-o, --output <OUTPUT_LOCATION>            
--compression <COMPRESSION>      [default: fast]         
--parents                    
--children                   
--no-compress                
--exclude                    
  -h, --help                       Print help
  -V, --version                    Print version
```

`--parents`: This will extract all the reads classified at all taxons between the root and the specified `--taxid`

`--children`: This will extract all the reads classified as decendents or subtaxa of `--taxid` (Including the taxid)

`--compression`: This defines the compression mode of the output `fastq.gz` file - fast / default / best

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

- 0.2.2

## Changelog

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
