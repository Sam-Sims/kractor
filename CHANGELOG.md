# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] - 2025-08-12

### Added
- Added a `reads_extracted_per_taxon` field to to summary report (#28)
- Added a `proportion_extracted` field to summary report (#28)
- Added the version to summary report (#28)
- Added an output format (`fasta` or `fastq`) field to the summary report (#28)
- Added a `--verbose` flag (in addition to the existing `-v`)

### Changed
- Removed `-O` for compression type, now uses `--compression-format` for clarity.
- Removed `-l` for compression level, now uses `--compression-level` for clarity.
- Renamed `--json-report` to `--summary`
- Improved the JSON report format to make it easier to read by removing `Paired` and `Single` fields and instead having a simple `total_reads_in` and `total_reads_out` field.

### Fixed
- Removed duplicate log message for taxon IDs identified
- Clippy warnings

## [1.0.1] - 2025-06-28

### Fixed
- Create subdirectories specified in output path if they don't exist (#24)
- Add output validation to prevent overwriting existing files (#25)

## [1.0.0] - 2025-04-16

### Added
- Better error handling with color-eyre (#16)
- Proper panic handling (#16)
- Tests for most functions
- JSON report with accurate read count information (#15)

### Changed
- Keep fastq parsing in bytes instead of converting to String (#17)
- Optimized functions to take `&str` instead of `String` (#21)
- Migrated to crossbeam scoped channels and refactored threading code
- Refactored JSON output and removed lazy_static dependency

### Fixed
- Root node no longer added to tree twice (#22)

## [0.4.0] - 2023-10-06

### Added
- JSON report included in stdout upon successful completion (can be disabled with `--no-json`)

### Changed
- Project renamed

## [0.3.0] - 2023-09-22

### Added
- Support for paired-end files
- Output FASTA file with `--output-fasta` option

### Changed
- Major refactor to use Noodles for fastq parsing
- Switched to Niffler for compression handling
- Streamlined arguments related to compression types/level

### Improved
- Logging functionality

## [0.2.3] - 2023-09-02

### Changed
- Code optimizations

## [0.2.2]

### Added
- `--no-compress` flag to output standard plaintext fastq files
- `--exclude` option to exclude specified reads (works with `--children` and `--parents`)
- Internal documentation (docstrings)

### Improved
- Increased verbosity of user outputs

## [0.2.1]

### Fixed
- Reduced memory usage

## [0.2.0]

### Added
- Automatic detection and handling of gz and plain files
- `--compression` argument to select compression type
- `--children` and `--parents` options to save children and parents based on kraken report

### Changed
- Integrated zlib-ng for faster gzip handling

## [0.1.0]

### Added
- Initial release