use anyhow::Result;
use chrono::Local;
use clap::Parser;
use env_logger::{fmt::Color, Builder};
use log::{info, LevelFilter};
use serde::{Deserialize, Serialize};
use std::io::prelude::*;

pub mod extract;
pub mod parsers;
pub use crate::cli::Cli;
pub mod cli;

use extract::{process_paired_end, process_single_end};
use parsers::kraken::process_kraken_output;

#[derive(Serialize, Deserialize)]
struct Summary {
    taxon_count: usize,
    taxon_ids: Vec<i32>,
    reads_output: ReadCounts,
    input_format: String,
}

#[derive(Serialize, Deserialize)]
#[serde(untagged)]
enum ReadCounts {
    Single {
        total: usize,
    },
    Paired {
        total: usize,
        read1: usize,
        read2: usize,
    },
}

/// Initializes and configures the logger.
///
/// This function sets up the logger for the application When verbosity is enabled, log messages
/// at the `Debug` level and above are output, else "Info" level is output.
///
/// # Arguments
///
/// * `verbose` - A boolean flag indicating whether verbose logging should be enabled.
fn logger(verbose: bool) {
    let level_filter = if verbose {
        LevelFilter::Debug
    } else {
        LevelFilter::Info
    };

    Builder::new()
        .format(|buf, record| {
            let mut style = buf.style();
            style.set_color(match record.level() {
                log::Level::Trace => Color::Magenta,
                log::Level::Debug => Color::Blue,
                log::Level::Info => Color::Green,
                log::Level::Warn => Color::Yellow,
                log::Level::Error => Color::Red,
            });

            writeln!(
                buf,
                "{} [{}] - {}",
                Local::now().format("[%H:%M:%S]"),
                style.value(record.level()),
                record.args()
            )
        })
        .filter(None, level_filter)
        .init();
}

fn main() -> Result<()> {
    let args = Cli::parse();

    //init logging
    logger(args.verbose);

    //Validate input/output args
    args.validate_input()?;

    //collect the taxon IDs to save and map those to read IDs
    let taxon_ids_to_save =
        extract::collect_taxons_to_save(&args.report, args.children, args.parents, args.taxid)?;
    let reads_to_save = process_kraken_output(&args.kraken, args.exclude, &taxon_ids_to_save)?;

    //check if paired-end reads are provided
    let paired = args.input.len() == 2;
    let input_format = if paired { "paired" } else { "single" };
    let read_counts = if paired {
        info!("Detected two input files. Assuming paired-end reads");
        let (total_reads_output_pair1, total_reads_output_pair2) = process_paired_end(
            &reads_to_save,
            args.input,
            args.output,
            args.output_type,
            args.compression_level,
            args.output_fasta,
        )?;

        ReadCounts::Paired {
            total: total_reads_output_pair1 + total_reads_output_pair2,
            read1: total_reads_output_pair1,
            read2: total_reads_output_pair2,
        }
    } else {
        info!("Detected one input file. Assuming single-end reads");
        let total_reads_output = process_single_end(
            &reads_to_save,
            args.input,
            args.output,
            args.output_type,
            args.compression_level,
            args.output_fasta,
        )?;

        ReadCounts::Single {
            total: total_reads_output,
        }
    };

    let summary = Summary {
        taxon_count: taxon_ids_to_save.len(),
        taxon_ids: taxon_ids_to_save,
        reads_output: read_counts,
        input_format: input_format.to_string(),
    };

    info!("Complete!");

    if !args.no_json {
        let json = serde_json::to_string_pretty(&summary)?;
        println!("{}", json);
    }

    Ok(())
}
