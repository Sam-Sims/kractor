use crate::parsers::fastx::{parse_fastq, write_output_fasta, write_output_fastq};
use crate::parsers::kraken::{build_tree_from_kraken_report, extract_children, extract_parents};
use anyhow::{anyhow, bail, Context, Result};
use crossbeam::{channel, thread};
use fxhash::FxHashSet;
use log::{debug, info};
use noodles::fastq;
use std::path::PathBuf;

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
pub fn process_single_end(
    reads_to_save: &FxHashSet<Vec<u8>>,
    input: Vec<PathBuf>,
    output: Vec<PathBuf>,
    compression_type: Option<niffler::Format>,
    compression_level: niffler::Level,
    fasta: bool,
) -> Result<usize> {
    thread::scope(|scope| -> Result<usize> {
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
        let total_reads_output = writer
            .join()
            .map_err(|_| anyhow!("Writer thread panicked"))??;
        Ok(total_reads_output)
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
pub fn process_paired_end(
    reads_to_save: &FxHashSet<Vec<u8>>,
    input: Vec<PathBuf>,
    output: Vec<PathBuf>,
    compression_type: Option<niffler::Format>,
    compression_level: niffler::Level,
    fasta: bool,
) -> Result<(usize, usize)> {
    thread::scope(|scope| -> Result<(usize, usize)> {
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
        let total_reads_output_pair1 = writer1
            .join()
            .map_err(|_| anyhow!("Writer thread for file1 panicked"))??;
        let total_reads_output_pair2 = writer2
            .join()
            .map_err(|_| anyhow!("Writer thread for file2 panicked"))??;
        Ok((total_reads_output_pair1, total_reads_output_pair2))
    })
    .map_err(|_| anyhow!("Thread communication error",))?
}

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
pub fn collect_taxons_to_save(
    report: &Option<PathBuf>,
    children: bool,
    parents: bool,
    taxids: Vec<i32>,
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
        let (nodes, taxon_map) = build_tree_from_kraken_report(&taxids, report_path)?;
        if children {
            debug!("Extracting children");
            let mut children = Vec::new();
            for taxid in &taxids {
                if let Some(&node_index) = taxon_map.get(taxid) {
                    extract_children(&nodes, node_index, &mut children)?;
                } else {
                    return Err(anyhow::anyhow!(
                        "Taxon ID {} not found in taxonomy map",
                        taxid
                    ));
                }
            }
            taxon_ids_to_save.extend(&children);
        } else if parents {
            debug!("Extracting parents");
            for taxid in &taxids {
                taxon_ids_to_save.extend(extract_parents(&taxon_map, &nodes, *taxid)?);
            }
        } else {
            taxon_ids_to_save.extend(&taxids);
        }
    } else {
        debug!(
            "No kraken report provided - extracting reads for taxon ID {:?} only",
            taxids
        );
        taxon_ids_to_save.extend(&taxids);
    }

    taxon_ids_to_save.sort_unstable();
    taxon_ids_to_save.dedup();

    debug!("Taxon IDs identified: {:?}", taxon_ids_to_save);
    if taxon_ids_to_save.is_empty() {
        bail!("No taxon IDs were identified for extraction");
    }
    Ok(taxon_ids_to_save)
}
