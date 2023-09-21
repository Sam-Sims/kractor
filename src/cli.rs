use clap::Parser;
use log::{debug, error, info, trace, warn};

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
    #[arg(short = 'k', long = "kraken", required = true)]
    pub kraken: String,
    // Kraken2 report file
    #[arg(short = 'r', long = "report")]
    pub report: Option<String>,
    // Taxid to extract reads for
    #[arg(short = 't', long = "taxid", required = true)]
    pub taxid: i32,
    // Compression type
    #[arg(short = 'O', long = "output-type", value_parser(validate_compression))]
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
    // #[arg(long)]
    // pub output_fasta: bool,
    // Verbose
    #[arg(short)]
    pub verbose: bool,
}

impl Cli {
    pub fn validate_input(&self) {
        let in_count = self.input.len();
        let out_count = self.output.len();
        if in_count > 2 {
            error!("Too many input files specified. Only 1 or 2 are allowed.");
            std::process::exit(1);
        }
        if out_count > 2 {
            error!("Too many output files specified. Only 1 or 2 are allowed.");
            std::process::exit(1);
        }
        if in_count == 2 && out_count == 1 {
            error!("Two input files specified but only one output file specified.");
            std::process::exit(1);
        }
        if in_count == 1 && out_count == 2 {
            error!("One input file specified but two output files specified.");
            std::process::exit(1);
        }
    }
}

fn validate_compression(s: &str) -> Result<niffler::compression::Format, String> {
    match s {
        "g" => Ok(niffler::compression::Format::Gzip),
        "b" => Ok(niffler::compression::Format::Bzip),
        "l" => Ok(niffler::compression::Format::Lzma),
        "z" => Ok(niffler::compression::Format::Zstd),
        "u" => Ok(niffler::compression::Format::No),
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
