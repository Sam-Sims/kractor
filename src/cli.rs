use clap::Parser;

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
    #[arg(long)]
    #[arg(long)]
    pub output: String,
    #[arg(long)]
    pub output2: Option<String>,
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
