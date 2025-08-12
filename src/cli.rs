use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(
    version,
    about = "Extract reads from a FASTQ file based on taxonomic classification via Kraken2.",
    author = "Sam Sims"
)]
pub struct Cli {
    /// Input file path(s). Accepts up to 2 files (for paired-end reads).
    #[arg(short = 'i', long = "input", num_args(0..=2), required = true)]
    pub input: Vec<PathBuf>,
    /// Output file path(s). Accepts up to 2 files (for paired-end reads).
    #[arg(short = 'o', long = "output", num_args(0..=2), required = true)]
    pub output: Vec<PathBuf>,
    /// Kraken2 stdout file path.
    #[arg(short = 'k', long = "kraken", required = true)]
    pub kraken: PathBuf,
    /// Kraken2 report file path.
    #[arg(short = 'r', long = "report", required_if_eq_any([("parents", "true"), ("children", "true")]))]
    pub report: Option<PathBuf>,
    /// One or more taxon IDs to extract reads for.
    #[arg(short = 't', long = "taxid", required = true, num_args(1..))]
    pub taxid: Vec<i32>,
    /// Include all parent taxon IDs in the output. Requires a Kraken2 report file.
    #[arg(short = 'p', long, action)]
    pub parents: bool,
    /// Include all child taxon IDs in the output. Requires a Kraken2 report file.
    #[arg(short = 'c', long, action)]
    pub children: bool,
    /// Compression format for output files (gz, bz2). Overides the inferred format.
    #[arg(long = "compression-format", value_parser(validate_compression))]
    pub output_type: Option<niffler::compression::Format>,
    /// Compression level (1-9).
    #[arg(
        long = "compression-level",
        default_value = "2",
        value_parser(validate_compression_level)
    )]
    pub compression_level: niffler::Level,
    /// Exclude specified taxon IDs from the output.
    #[arg(long)]
    pub exclude: bool,
    /// Output results in FASTA format
    #[arg(long, action)]
    pub output_fasta: bool,
    /// Enable a JSON summary output written to stdout.
    #[arg(long = "summary")]
    pub summary: bool,
    /// Enable verbose output.
    #[arg(short, long)]
    pub verbose: bool,
}

fn validate_compression(s: &str) -> Result<niffler::compression::Format, String> {
    match s {
        "gz" => Ok(niffler::compression::Format::Gzip),
        "bz2" => Ok(niffler::compression::Format::Bzip),
        "none" => Ok(niffler::compression::Format::No),
        _ => Err(format!("Unknown compression type: {s}")),
    }
}

fn validate_compression_level(s: &str) -> Result<niffler::Level, String> {
    match s.parse::<u32>() {
        Ok(1) => Ok(niffler::Level::One),
        Ok(2) => Ok(niffler::Level::Two),
        Ok(3) => Ok(niffler::Level::Three),
        Ok(4) => Ok(niffler::Level::Four),
        Ok(5) => Ok(niffler::Level::Five),
        Ok(6) => Ok(niffler::Level::Six),
        Ok(7) => Ok(niffler::Level::Seven),
        Ok(8) => Ok(niffler::Level::Eight),
        Ok(9) => Ok(niffler::Level::Nine),
        _ => Err(format!(
            "Unknown compression level: {s} Try a value between 1-9"
        )),
    }
}
