![GitHub release (with filter)](https://img.shields.io/github/v/release/sam-sims/Kractor)
![crates.io](https://img.shields.io/crates/v/kractor
)
![](https://anaconda.org/bioconda/kractor/badges/version.svg)
[![test](https://github.com/Sam-Sims/kractor/actions/workflows/test.yaml/badge.svg?branch=main)](https://github.com/Sam-Sims/kractor/actions/workflows/test.yaml)
[![check](https://github.com/Sam-Sims/kractor/actions/workflows/check.yaml/badge.svg?branch=main)](https://github.com/Sam-Sims/kractor/actions/workflows/check.yaml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15761837.svg)](https://doi.org/10.5281/zenodo.15761837)

# kractor

**kra**ken extr**actor**

Kractor extracts reads from FASTQ or FASTA files (optionally `.gz`/`.bz2` compressed) using taxonomic classifications
from Kraken2. It supports single- and paired-end reads, can optionally include taxonomic parents or children, and uses
minimal memory (~4.5 MB for a 17 GB FASTQ file).

The output is a FASTQ or FASTA file containing all reads classified as the specified taxa.

Kractor offers a significant speed-up over KrakenTools for both paired and unpaired reads.

Performance vs KrakenTools:

- Paired compressed FASTQ: ~21× faster
- Paired uncompressed FASTQ: ~10× faster
- Unpaired: ~4× faster (compressed or uncompressed)

For additional details, refer to the [benchmarks](benchmarks/benchmarks.md)

## Motivation

Provides similar functionality to the [KrakenTools](https://github.com/jenniferlu717/KrakenTools) `extract_kraken_reads`
python script.

However the main motivation was to enhance speed when processing multiple, large FASTQ files - and as a way to learn
Rust.

## Installation

### Binaries:

Precompiled binaries for Linux, MacOS and Windows are attached to the latest
release.

### Conda:

![](https://anaconda.org/bioconda/kractor/badges/platforms.svg)

```
conda install -c bioconda kractor
```

### Docker:

A docker image is available on [quay.io](https://quay.io/repository/biocontainers/kractor)

```
docker pull quay.io/biocontainers/kractor
```

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

```bash
Extract reads from a FASTQ or FASTA file based on taxonomic classification via Kraken2.

Usage: kractor [OPTIONS] --input [<INPUT>...] --output [<OUTPUT>...] --kraken <KRAKEN> --taxid <TAXID>...

Options:
  -i, --input [<INPUT>...]
          Input file path(s). Accepts up to 2 files (for paired-end reads)
  -o, --output [<OUTPUT>...]
          Output file path(s). Accepts up to 2 files (for paired-end reads)
  -k, --kraken <KRAKEN>
          Kraken2 stdout file path
  -r, --report <REPORT>
          Kraken2 report file path
  -t, --taxid <TAXID>...
          One or more taxon IDs to extract reads for
  -p, --parents
          Include all parent taxon IDs in the output. Requires a Kraken2 report file
  -c, --children
          Include all child taxon IDs in the output. Requires a Kraken2 report file
      --compression-format <OUTPUT_TYPE>
          Compression format for output files (gz, bz2). Overrides the inferred format
      --compression-level <COMPRESSION_LEVEL>
          Compression level (1-9) [default: 2]
      --exclude
          Exclude specified taxon IDs from the output
      --output-format <OUTPUT_FORMAT>
          Output sequence format (auto, fasta, fastq) [default: auto] [possible values: auto, fasta, fastq]
      --summary
          Enable a JSON summary output written to stdout
      --no-header-detect
          Disable detection and skipping of any header lines in the Kraken2 report
  -v, --verbose
          Enable verbose output
  -h, --help
          Print help
  -V, --version
          Print version
  ```

### Examples:

```bash
# Extract reads classified as E. coli from single end reads
kractor -i sample.fastq -o extracted.fastq -k kraken_output.txt -t 562

# Extract from paired end reads
kractor -i sample_R1.fastq -i sample_R2.fastq -o extracted_R1.fastq -o extracted_R2.fastq -k kraken_output.txt -t 562

# Extract multiple taxids (Bacillaceae and Listeriaceae)
kractor -i sample.fastq -o extracted.fastq -k kraken_output.txt -t 186817 186820

# Extract all children of Enterobacteriaceae family (requires kraken report)
kractor -i sample.fastq -o extracted.fastq -k kraken_output.txt -r kraken_report.txt -t 543 --children

# Extract everything EXCEPT viral reads (using --exclude)
kractor -i sample.fastq -o extracted.fastq -k kraken_output.txt -t 10239 --exclude

# Output FASTA format instead of FASTQ
kractor -i sample.fastq -o extracted.fasta -k kraken_output.txt -t 562 --output-format fasta
```

### Summary statistics

Use `--summary` to get summary statistics (output to stdout on completion)

```json
{
  "kractor_version": "5.0.0",
  "input_layout": "single",
  "input_sequence_format": "fastq",
  "output_sequence_format": "fasta",
  "requested_taxon_ids": [
    0,
    1,
    999999999
  ],
  "matched_taxon_ids": [
    0,
    1
  ],
  "requested_taxon_ids_not_found": [
    999999999
  ],
  "total_input_records": 3491078,
  "total_output_records": 747237,
  "extraction_fraction": 0.2140419091180432,
  "assigned_reads_per_taxon": {
    "0": 745591,
    "1": 1646
  }
}
```

Fields:

- `kractor_version`: Version of kractor that produced the summary.
- `input_layout`: `single` or `paired` input mode.
- `input_sequence_format`: Input sequence format, `fastq` or `fasta`.
- `output_sequence_format`: Output sequence format, `fastq` or `fasta`.
- `requested_taxon_ids`: Taxon IDs requested with `--taxid`.
- `matched_taxon_ids`: Requested taxon IDs found in the Kraken report/output (includes parents/children if used).
- `requested_taxon_ids_not_found`: Requested taxon IDs that were not found in the Kraken report/output.
- `total_input_records`: Total records parsed from the input file(s).
- `total_output_records`: Total records written to the output file(s).
- `extraction_fraction`: `total_output_records / total_input_records`.
- `assigned_reads_per_taxon`: Number of directly assigned reads per matched taxon ID (0 indicates no direct assignments, but present due to children/parents).

### Arguments:

### Required:

#### Input

`-i, --input`

Specifies one input FASTA or FASTQ file, or two input FASTQ files for paired-end extraction. Files may be uncompressed or compressed (`gz`, `bz2`).
Paired-end reads can be specified by:

Using `--input` twice: `-i <R1_fastq_file> -i <R2_fastq_file>`

Using `--input` once but passing both files: `-i <R1_fastq_file> <R2_fastq_file>`

#### Output

`-o, --output`

Specifies the output file(s) for extracted reads, matching the order of the input files.
Compression type is inferred from the file extension (`.gz`, `.bz2`). If not recognised, output will be uncompressed. Sequence format defaults to the input format unless `--output-format` is set.

#### Kraken Output

`-k, --kraken`

Path to
the [Standard Kraken Output Format file](https://github.com/DerrickWood/kraken2/wiki/Manual#standard-kraken-output-format),
containing taxonomic classification of read IDs.

#### Taxid

`-t, --taxid`

One or more taxonomic IDs to extract.

For example: `-t 1 2 10`

Each taxonomic id is affected by `--exclude`, `--parents`, and `--children` if those options are used.

Taxonomic ids do not need to be present in a given report. This may be useful when running kractor in a wrapper script
for several FASTQ/FASTA files and just want to extract a set of taxon ids from them all - without caring if they are present
or not.

### Optional:

#### Compression format

`--compression-format`

Manually set output compression format, overriding what is inferred from file names.

Valid values:

- `gz` – gzip compression
- `bz2` – bzip2 compression
- `none` – no compression

#### Compression level

`--compression-level`

Set compression level (1–9).

- 1 = fastest, largest file
- 9 = slowest, smallest file

Default: 2 (balance of speed and size)

#### Output sequence format

`--output-format`

Set the output sequence format.

Valid values:

- `auto` – use the input sequence format (default)
- `fasta` – write FASTA output
- `fastq` – write FASTQ output

#### Kraken Report

`-r, --report`

Path to the [Kraken2 report file](https://github.com/DerrickWood/kraken2/wiki/Manual#sample-report-output-format).
Required if using `--parents` or `--children`.

The first line is automatically treated as a header if it looks non-numeric; use `--no-header-detect` to force parsing
from the very first line. Parsing errors will include the report line number and offending line to help spot format
issues.

#### Parents

`--parents`

Include reads classified between the root and the specified `--taxid`. Requires `--report`.

#### Children

`--children`

Include reads classified at the given taxid and all its descendant taxa. Requires `--report`.

#### Exclude

`--exclude`

Extract all reads except those matching the given taxids. Can be combined with `--parents` or `--children`.

#### JSON report

`--summary`

Write a JSON report to stdout after processing.

## Citation

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15761837.svg)](https://doi.org/10.5281/zenodo.15761837)

```
Sam Sims. (2025). Sam-Sims/kractor. Zenodo. https://doi.org/10.5281/zenodo.15761837
```
