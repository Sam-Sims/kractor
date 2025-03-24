use anyhow::{anyhow, Result};
use clap::Parser;

#[derive(Parser, Debug)]
#[command(
    version,
    about = "Extract reads from a FASTQ file based on taxonomic classification via Kraken2."
)]
pub struct Cli {
    // Fastq file(s)
    #[arg(short = 'i', long = "input", num_args(0..=2), required = true, value_parser(check_input_exists))]
    pub input: Vec<String>,
    // Output file(s)
    #[arg(short = 'o', long = "output", num_args(0..=2), required = true)]
    pub output: Vec<String>,
    // Kraken2 output file
    #[arg(
        short = 'k',
        long = "kraken",
        required = true,
        value_parser(check_input_exists)
    )]
    pub kraken: String,
    // Kraken2 report file
    #[arg(short = 'r', long = "report", value_parser(check_input_exists), required_if_eq_any([("parents", "true"), ("children", "true")]))]
    pub report: Option<String>,
    // Taxid to extract reads for
    #[arg(short = 't', long = "taxid", required = true)]
    pub taxid: i32,
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

impl Cli {
    pub fn validate_input(&self) -> Result<()> {
        let in_count = self.input.len();
        let out_count = self.output.len();
        if in_count > 2 {
            return Err(anyhow!(
                "Too many input files specified. Only 1 or 2 are allowed."
            ));
        }
        if out_count > 2 {
            return Err(anyhow!(
                "Too many output files specified. Only 1 or 2 are allowed."
            ));
        }
        if in_count == 2 && out_count == 1 {
            return Err(anyhow!(
                "Two input files specified but only one output file specified."
            ));
        }
        if in_count == 1 && out_count == 2 {
            return Err(anyhow!(
                "One input file specified but two output files specified."
            ));
        }
        Ok(())
    }
}

fn validate_compression(s: &str) -> Result<niffler::compression::Format, String> {
    match s {
        "gz" => Ok(niffler::compression::Format::Gzip),
        "bz2" => Ok(niffler::compression::Format::Bzip),
        "none" => Ok(niffler::compression::Format::No),
        _ => Err(format!("Unknown compression type: {}", s)),
    }
}

fn check_input_exists(s: &str) -> Result<String, String> {
    if std::path::Path::new(s).exists() {
        Ok(s.to_string())
    } else {
        Err(format!("File does not exist: {}", s))
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
