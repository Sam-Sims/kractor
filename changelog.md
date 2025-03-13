
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
