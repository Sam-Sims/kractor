use anyhow::{anyhow, Result};
use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(
    version,
    about = "Extract reads from a FASTQ file based on taxonomic classification via Kraken2.",
    author = "Sam Sims"
)]
pub struct Cli {
    // Fastq file(s)
    #[arg(short = 'i', long = "input", num_args(0..=2), required = true)]
    pub input: Vec<PathBuf>,
    // Output file(s)
    #[arg(short = 'o', long = "output", num_args(0..=2), required = true)]
    pub output: Vec<PathBuf>,
    // Kraken2 output file
    #[arg(short = 'k', long = "kraken", required = true)]
    pub kraken: PathBuf,
    // Kraken2 report file
    #[arg(short = 'r', long = "report", required_if_eq_any([("parents", "true"), ("children", "true")]))]
    pub report: Option<PathBuf>,
    // Taxid to extract reads for
    #[arg(short = 't', long = "taxid", required = true, num_args(1..))]
    pub taxid: Vec<i32>,
    // Compression type
    #[arg(
        short = 'O',
        long = "compression-type",
        value_parser(validate_compression)
    )]
    pub output_type: Option<niffler::compression::Format>,
    //Compression level
    #[arg(
        short = 'l',
        long = "level",
        default_value = "2",
        value_parser(validate_compression_level)
    )]
    pub compression_level: niffler::Level,
    // Extract reads from parents
    #[arg(long, action)]
    pub parents: bool,
    // Extract reads from children
    #[arg(long, action)]
    pub children: bool,
    // Exclude reads matching taxid
    #[arg(long)]
    pub exclude: bool,
    // Output reads in FASTA format
    #[arg(long, action)]
    pub output_fasta: bool,
    // Dont output json
    #[arg(long = "no-json")]
    pub no_json: bool,
    // Verbose
    #[arg(short)]
    pub verbose: bool,
}

fn validate_compression(s: &str) -> Result<niffler::compression::Format, String> {
    match s {
        "gz" => Ok(niffler::compression::Format::Gzip),
        "bz2" => Ok(niffler::compression::Format::Bzip),
        "none" => Ok(niffler::compression::Format::No),
        _ => Err(format!("Unknown compression type: {}", s)),
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
            "Unknown compression level: {} Try a value between 1-9",
            s
        )),
    }
}
