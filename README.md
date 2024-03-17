[![Release](https://github.com/Sam-Sims/Kractor/actions/workflows/release.yaml/badge.svg)](https://github.com/Sam-Sims/Kractor/actions/workflows/release.yaml)
![GitHub release (with filter)](https://img.shields.io/github/v/release/sam-sims/Kractor)
![crates.io](https://img.shields.io/crates/v/kractor
)
![Docker Pulls](https://img.shields.io/docker/pulls/samsims/kractor)


# Kractor

**kra**ken extr**actor**

Kractor extracts sequencing reads based on taxonomic classifications obtained via [Kraken2](https://github.com/DerrickWood/kraken2). It consumes paired or unpaired `fastq[.gz/.bz]` files as input alongisde a Kraken2 standard output. It can optionally consume a Kraken2 report to extract all taxonomic parents and children of a given taxid. Fast by default, it outputs `fast[q/a]` files, that can optionally be compressed.

Kractor significantly enhances processing speed compared to KrakenTools for both paired and unpaired reads. Paired reads are processed approximately 21x quicker for compressed fastqs and 10x quicker for uncompressed. Unpaired reads are approximately 4x faster for both compressed and uncompressed inputs.

 For additional details, refer to the [benchmarks](benchmarks/benchmarks.md)

## Motivation

Heavily inspired by the great [KrakenTools](https://github.com/jenniferlu717/KrakenTools). 

At the time of writing KrakenTools operates as a single-threaded Python implementation which poses limitations in speed when processing large, paired-end fastq files. The main motivation was to enchance speed when parsing and extracting (writing) a large volume of reads - and also to learn rust!

## Installation

### Binaries:

Precompiled binaries for Linux, MacOS and Windows are attached to the latest release [0.4.0](https://github.com/Sam-Sims/Kractor/releases/tag/v0.4.0)

### Docker:

A docker image is available on [Docker Hub](https://hub.docker.com/r/samsims/kractor)

```bash
docker pull samsims/kractor
docker run samsims/kractor --help
```

Use `-v` to mount your input and output directories. A typical command might look like:

```bash
docker run -v /path/to/input:/input -v /path/to/output:/output samsims/kractor -k /input/<kraken_output> -i /input/<fastq_file> -t <taxonomic_id> -o /output/<output_fastq>
````

### Cargo:
Requires [cargo](https://www.rust-lang.org/tools/install)
```
cargo install kractor
```

### Build from source:

#### Install rust toolchain:

To install please refer to the rust documentation: [docs](https://www.rust-lang.org/tools/install)

#### Clone the repository:

```bash
git clone https://github.com/Sam-Sims/Kractor
```

#### Build and add to path:

```bash
cd Kractor
cargo build --release
export PATH=$PATH:$(pwd)/target/release
```

All executables will be in the directory Kractor/target/release.

## Usage
![Alt text](screenshot.png)
### Basic Usage:

```bash
kractor -k <kraken_output> -i <fastq_file> -t <taxonomic_id> -o <output_file> > kractor_report.json
```
Or, if you have paired-end illumina reads:
```bash
kractor -k <kraken_output> -i <R1_fastq_file> -i <R2_fastq_file> -t <taxonomic_id> -o <R1_output_file> -o <R2_output_file>
```
If you want to extract all children of a taxon:
```bash
kractor -k <kraken_output> -r <kraken_report> -i <fastq_file> -t <taxonomic_id> --children -o <output_file>
```

### Arguments:

### Required:

#### Input

`-i, --input`

This option will specify the input files containing the reads you want to extract from. They can be compressed - (`gz`, `bz2`). Paired end reads can be specified by:

Using `--input` twice: `-i <R1_fastq_file> -i <R2_fastq_file>`

Using `--input` once but passing both files: `-i <R1_fastq_file> <R2_fastq_file>`

  This means that bash wildcard expansion works: `-i *.fastq`

#### Output

`-o, --output`

This option will specify the output files containing the extracted reads. The order of the output files is assumed to be the same as the input. 

By default the compression will be inferred from the output file extension for supported file types (`gz`, `bz`). If the output type cannot be inferred, plaintext will be output.

#### Kraken Output

`-k, --kraken`

This option will specify the path to the Kraken2 output containing taxonomic classification of read IDs.

#### Taxid

`-t, --taxid`

This option will specify the taxon ID for reads you want to extract.

### Optional:

#### Output type

`-O, --output-type`

This option will manually set the compression mode used for the output file and will override the type inferred from the output path. 

Valid values are:

- `gz` to output gz
- `bz2` to output bz2
- `none` to not apply compresison

#### Compression level

`-l, --level`

This option will set the compression level to use if compressing the output. Should be a value between 1-9 with 1 being the fastest but largest file size and 9 is for slowest, but best file size. By default this is set at 2 as it is a good trade off for speed/filesize.

#### Output fasta

`--output-fasta`

This option will output a fasta file, with read ids as headers.

#### Kraken Report

`-r, --report`

This option specifies the path to the report file generated by Kraken2. If you want to use `--parents` or `--children` then is argument is required.

#### Parents

`--parents`

This will extract reads classified at all taxons between the root and the specified `--taxid`.

#### Children

`--children`

This will extract all the reads classified as decendents or subtaxa of `--taxid` (Including the taxid).

#### Exclude

`--exclude`

This will output every read except those matching the taxid. Works with `--parents` and `--children`

#### Skip report

`--no-json`

This will skip the json report that is output to stdout upon programme completion.

## Future plans

- [x] Support unzipped fastq files
- [x] Support paired end FASTQ files
- [x] `--include-parents` and `--include-children` arguments
- [ ] Supply multiple taxonomic IDs to extract
- [x] Exclude taxonomic IDs
- [x] `--compression-mode`
- [x] More verbose output
- [x] Benchmarks
- [x] Output fasta
- [x] Output non `gz`
- [ ] Tests

## Version

- 0.4.0

## Changelog

### 0.4.0
- Json report including in stdout upon successful completion (can be disabled with --no-json)
- Renamed

### 0.3.0

- Support for paired-end files
- Major under the hood changes to use Noodles for fastq parsing, and Niffler to handle compression (RIP my own code)
- Output a fasta file with `--output-fasta`
- Streamline arguments related to compression types/level
- Improved logging

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
