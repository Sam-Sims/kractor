use clap::Parser;
use log::{debug, error, info, trace, warn};

#[derive(Parser, Debug)]
#[command(
    version,
    about = "Extract reads from a FASTQ file based on taxonomic classification via Kraken2."
)]
pub struct Cli {
    // Fastq file to extract reads from
    #[arg(short = 'i', long = "input", num_args(0..=2), required = true)]
    pub input: Vec<String>,
    #[arg(short, long)]
    pub kraken: String,
    #[arg(short, long)]
    pub taxid: i32,
    #[arg(short, long)]
    pub report: Option<String>,
    #[arg(short = 'o', long = "output", num_args(0..=2), required = true)]
    pub output: Vec<String>,
    #[arg(long, default_value = "fast")]
    pub compression_mode: Option<String>,
    #[arg(long, action)]
    pub parents: bool,
    #[arg(long, action)]
    pub children: bool,
    #[arg(long)]
    pub no_compress: bool,
    #[arg(long)]
    pub exclude: bool,
    #[arg(long)]
    pub output_fasta: bool,
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
