pub use crate::cli::Cli;
use anyhow::{anyhow, bail, Context, Result};
use chrono::Local;
use clap::Parser;
use crossbeam::channel;
use crossbeam::thread;
use env_logger::{fmt::Color, Builder};
use fxhash::FxHashSet;
use log::{debug, info, trace, LevelFilter};
use noodles::fastq;
use std::io::prelude::*;
use std::path::PathBuf;

pub mod parsers;

mod cli;
use parsers::fastx::{parse_fastq, write_output_fasta, write_output_fastq};
use parsers::kraken::{
    build_tree_from_kraken_report, extract_children, extract_parents, process_kraken_output,
};

/// Collects taxon IDs to save.
///
/// This function determines what taxon IDs need to be saved from the kraken output.
/// If a Kraken report is specified, it builds a tree of all taxons in the report and extracts taxon IDs based
/// on if --children or --parent are supplied. If no report is provided, the function returns only the given taxon ID
/// in the list of taxon IDs to save.
///
/// # Arguments
///
/// * `args` - The Args structure containing command-line arguments.
///
/// # Returns
///
/// A vector of taxon IDs that need to be saved.
fn collect_taxons_to_save(
    report: &Option<PathBuf>,
    children: bool,
    parents: bool,
    taxid: i32,
) -> Result<Vec<i32>> {
    let mut taxon_ids_to_save = Vec::new();

    // I dont think we will reach this code ever since clap should catch this - but in case it doesnt
    if (parents || children) && report.is_none() {
        return Err(anyhow::anyhow!(
            "Report required when parents or children is enabled"
        ));
    }

    if let Some(report_path) = report {
        info!("Processing kraken report...");
        let (nodes, taxon_map) = build_tree_from_kraken_report(taxid, &report_path)?;
        if children {
            debug!("Extracting children");
            let mut children = Vec::new();
            extract_children(&nodes, taxon_map[&taxid], &mut children)?;
            taxon_ids_to_save.extend(&children);
        } else if parents {
            debug!("Extracting parents");
            taxon_ids_to_save.extend(extract_parents(&taxon_map, &nodes, taxid)?);
        } else {
            taxon_ids_to_save.push(taxid);
        }
    } else {
        debug!(
            "No kraken report provided - extracting reads for taxon ID {} only",
            taxid
        );
        taxon_ids_to_save.push(taxid);
    }
    debug!("Taxon IDs identified: {:?}", taxon_ids_to_save);
    if taxon_ids_to_save.is_empty() {
        bail!("No taxon IDs were identified for extraction");
    }
    Ok(taxon_ids_to_save)
}

/// Process single-end reads from a FASTQ file.
///
/// This function runs the processing of single-end reads from a FASTQ file.
/// It spawns two threads: one for reading and parsing the FASTQ file and another for writing
/// the selected reads to an output file.
///
/// # Arguments
///
/// * `reads_to_save` - A HashMap containing read IDs and their associated taxon IDs.
/// * `input` - A vector containing the paths to the input file.
/// * `output` - A vector containing the paths to the output file.
/// * `output_type` - The compression type to use for the output file.
/// * `compression_level` - The compression level to use for the output file.
/// * `fasta` - A boolean indicating whether the output should be in FASTA format.
fn process_single_end(
    reads_to_save: &FxHashSet<Vec<u8>>,
    input: Vec<PathBuf>,
    output: Vec<PathBuf>,
    compression_type: Option<niffler::Format>,
    compression_level: niffler::Level,
    fasta: bool,
) -> Result<()> {
    thread::scope(|scope| -> Result<()> {
        let (tx, rx) = channel::unbounded::<fastq::Record>();

        let reader = scope.spawn(|_| {
            let result = parse_fastq(&input[0], reads_to_save, &tx);
            drop(tx);
            result.with_context(|| format!("Failed to parse input file: {:?}", input[0]))
        });

        let writer = scope.spawn(|_| {
            if !fasta {
                write_output_fastq(rx, &output[0], compression_type, compression_level)
                    .with_context(|| format!("Failed to write output file: {:?}", output[0]))
            } else {
                write_output_fasta(rx, &output[0])
                    .with_context(|| format!("Failed to write output file: {:?}", output[0]))
            }
        });

        reader
            .join()
            .map_err(|_| anyhow!("Reader thread panicked"))??;
        writer
            .join()
            .map_err(|_| anyhow!("Writer thread panicked"))??;
        Ok(())
    })
    .map_err(|_| anyhow!("Thread communication error"))?
}

/// Process paired-end reads from FASTQ files.
///
/// This function runs the processing of paired-end reads from two FASTQ input files.
/// It spawns two threads for reading and parsing each input file and two threads for writing
/// selected reads to their respective output files.
///
/// # Arguments
///
/// * `reads_to_save` - A HashMap containing read IDs and their associated taxon IDs.
/// * `input` - A vector containing the paths to the two input files.
/// * `output` - A vector containing the paths to the two output files.
/// * `compression_type` - The compression type to use for the output files.
/// * `compression_level` - The compression level to use for the output files.
/// * `fasta` - A boolean indicating whether to output in FASTA format.
fn process_paired_end(
    reads_to_save: &FxHashSet<Vec<u8>>,
    input: Vec<PathBuf>,
    output: Vec<PathBuf>,
    compression_type: Option<niffler::Format>,
    compression_level: niffler::Level,
    fasta: bool,
) -> Result<()> {
    thread::scope(|scope| -> Result<()> {
        let (tx1, rx1) = channel::unbounded::<fastq::Record>();
        let (tx2, rx2) = channel::unbounded::<fastq::Record>();

        let reader1 = scope.spawn(|_| {
            let result = parse_fastq(&input[0], reads_to_save, &tx1);
            drop(tx1);
            result.with_context(|| format!("Failed to parse first input file: {:?}", input[0]))
        });

        let reader2 = scope.spawn(|_| {
            let result = parse_fastq(&input[1], reads_to_save, &tx2);
            drop(tx2);
            result.with_context(|| format!("Failed to parse second input file: {:?}", input[1]))
        });

        let writer1 = scope.spawn(|_| {
            if !fasta {
                write_output_fastq(rx1, &output[0], compression_type, compression_level)
                    .with_context(|| {
                        format!(
                            "Failed to write FASTQ output to first file: {:?}",
                            output[0]
                        )
                    })
            } else {
                write_output_fasta(rx1, &output[0]).with_context(|| {
                    format!(
                        "Failed to write FASTA output to first file: {:?}",
                        output[0]
                    )
                })
            }
        });

        let writer2 = scope.spawn(|_| {
            if !fasta {
                write_output_fastq(rx2, &output[1], compression_type, compression_level)
                    .with_context(|| {
                        format!(
                            "Failed to write FASTQ output to second file: {:?}",
                            output[1]
                        )
                    })
            } else {
                write_output_fasta(rx2, &output[1]).with_context(|| {
                    format!(
                        "Failed to write FASTA output to second file: {:?}",
                        output[1]
                    )
                })
            }
        });

        reader1
            .join()
            .map_err(|_| anyhow!("Reader thread for file1 panicked"))??;
        reader2
            .join()
            .map_err(|_| anyhow!("Reader thread for file2 panicked"))??;
        writer1
            .join()
            .map_err(|_| anyhow!("Writer thread for file1 panicked"))??;
        writer2
            .join()
            .map_err(|_| anyhow!("Writer thread for file2 panicked"))??;
        Ok(())
    })
    .map_err(|_| anyhow!("Thread communication error",))?
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
        collect_taxons_to_save(&args.report, args.children, args.parents, args.taxid)?;
    let reads_to_save = process_kraken_output(&args.kraken, args.exclude, &taxon_ids_to_save)?;

    //check if paired-end reads are provided
    let paired = args.input.len() == 2;
    if paired {
        info!("Detected two input files. Assuming paired-end reads");
        process_paired_end(
            &reads_to_save,
            args.input,
            args.output,
            args.output_type,
            args.compression_level,
            args.output_fasta,
        )?;
    } else {
        info!("Detected one input file. Assuming single-end reads");
        process_single_end(
            &reads_to_save,
            args.input,
            args.output,
            args.output_type,
            args.compression_level,
            args.output_fasta,
        )?;
    }

    info!("Complete!");

    if !args.no_json {
        // let stats = serde_json::json!({
        //     "taxon_id_count": *TAXON_ID_COUNT.lock()
        //         .map_err(|e| anyhow!("Failed to lock TAXON_ID_COUNT mutex: {}", e))?,
        //     "taxon_ids": *TAXON_IDS.lock()
        //         .map_err(|e| anyhow!("Failed to lock TAXON_IDS mutex: {}", e))?,
        //     "reads_in": *TOTAL_READS.lock()
        //         .map_err(|e| anyhow!("Failed to lock TOTAL_READS mutex: {}", e))?,
        //     "reads_out": *READS_TO_EXTRACT.lock()
        //         .map_err(|e| anyhow!("Failed to lock READS_TO_EXTRACT mutex: {}", e))?,
        //     "input_format": if paired { "paired-end" } else { "single-end" },
        //     "output_format": if args.output_fasta { "fasta" } else { "fastq" },
        // });

        // let stats_json =
        //     serde_json::to_string_pretty(&stats).context("Failed to serialize stats to JSON")?;
        //
        // println!("{}", stats_json);
    }

    Ok(())
}
