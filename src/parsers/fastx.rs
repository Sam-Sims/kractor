use anyhow::Context;
use anyhow::Result;
use crossbeam::channel::{Receiver, Sender};
use log::{debug, trace};
use noodles::fasta::record::{Definition, Sequence};
use noodles::{fasta, fastq};
use std::collections::HashSet;
use std::io::BufReader;
use std::path::Path;
use std::sync::Arc;
use std::time::{Duration, Instant};
use std::{fs, io};

/// Parse a FASTQ file and send reads to writer thread.
///
/// This function reads a fastq file and extracts read IDs and sequences.
/// It compares the read IDs against a given HashMap of read IDs (`reads_to_save`) and sends
/// the sequences of matching read IDs to the writer thread.
///
/// # Arguments
///
/// * `file_path` - A string containing the path to the FASTQ file.
/// * `reads_to_save` - A HashMap containing read IDs and their associated taxon IDs.
/// * `tx` - A Sender channel to send the parsed reads to the writer thread.
pub fn parse_fastq(
    file_path: &str,
    reads_to_save: Arc<HashSet<String>>,
    tx: &Sender<fastq::Record>,
) -> Result<()> {
    let mut num_reads = 0;

    let mut last_progress_update = Instant::now();
    const PROGRESS_UPDATE_INTERVAL: Duration = Duration::from_millis(1500);

    let (reader, format) = niffler::from_path(file_path)
        .with_context(|| format!("Failed to open fastq file: {}", file_path))?;
    debug!(
        "Detected input compression type for file {:?} as: {:?}",
        file_path, format
    );
    let reader = BufReader::new(reader);
    let mut fastq_reader = fastq::Reader::new(reader);

    for (record_idx, result) in fastq_reader.records().enumerate() {
        let record = result
            .with_context(|| format!("Error reading FASTQ record at position {}", record_idx))?;

        let read_name_bytes = record.name();
        let read_name_string = std::str::from_utf8(read_name_bytes).with_context(|| {
            format!("Invalid UTF-8 sequence in read name: {:?}", read_name_bytes)
        })?;

        if reads_to_save.contains(&read_name_string.to_string()) {
            tx.send(record)
                .context("Failed to send record to writer thread")?;
        }

        num_reads += 1;

        if last_progress_update.elapsed() >= PROGRESS_UPDATE_INTERVAL {
            trace!("Processed {} reads", num_reads);
            last_progress_update = Instant::now();
        }
    }

    Ok(())
}

/// Infer the compression format from a file path based on its extension.
///
/// This function takes a file path as input and checks its file extension to determine
/// the compression format. It supports all niffler compression  types
/// If the extension is not recognized, it defaults to no compression.
///
/// # Arguments
///
/// * `file_path` - A reference to a `String` containing the file path to analyze.
///
/// # Returns
///
/// Niffler compression format.
fn infer_compression(file_path: &str) -> niffler::compression::Format {
    let path = Path::new(&file_path);
    let ext = path.extension().unwrap().to_str().unwrap();
    match ext {
        "gz" => niffler::compression::Format::Gzip,
        "bz2" => niffler::compression::Format::Bzip,
        _ => niffler::compression::Format::No,
    }
}

/// Write fastq data to output file
///
/// This function takes received data from the provided Receiver channel corresponding to an inputfile and writes it to the specified output
/// file to an output file, optionally applying compression.
///
/// # Arguments
///
/// * `rx`: A Receiver from which parsed fastq lines will be received.
/// * `out_file`: A file representing the output file where data will be written.
/// * `output_type`: The compression type to use for the output file.
/// * `compression_level`: The compression level to use for the output file.
pub fn write_output_fastq(
    rx: Receiver<fastq::Record>,
    out_file: &str,
    output_type: Option<niffler::Format>,
    compression_level: niffler::Level,
) -> Result<()> {
    let compression_type = match output_type {
        Some(output_type) => {
            debug!("Output type overridden as: {:?}", output_type);
            output_type
        }
        None => {
            let inferred_type = infer_compression(out_file);
            debug!("Inferred output compression type as: {:?}", inferred_type);
            inferred_type
        }
    };

    debug!(
        "Output compression level specified as: {:?}",
        compression_level
    );
    debug!("Creating output file: {:?}", out_file);

    let out_file = fs::File::create(out_file)
        .with_context(|| format!("Failed to create output file: {}", out_file))?;

    let file_handle = Box::new(io::BufWriter::new(out_file));
    let writer = niffler::get_writer(file_handle, compression_type, compression_level)
        .context("Failed to create niffler writer")?;

    let mut fastq_writer = fastq::Writer::new(writer);

    for record in rx {
        fastq_writer
            .write_record(&record)
            .with_context(|| format!("Error writing FASTQ record: {:?}", record))?;
    }

    Ok(())
}

pub fn write_output_fasta(rx: Receiver<fastq::Record>, out_file: &str) -> Result<()> {
    debug!("Creating output file: {:?}", out_file);

    let out_file = fs::File::create(out_file)
        .with_context(|| format!("Failed to create output file: {}", out_file))?;

    let mut writer = fasta::Writer::new(out_file);

    for record in rx {
        let definition = Definition::new(
            std::str::from_utf8(record.name()).with_context(|| {
                format!("Invalid UTF-8 sequence in read name: {:?}", record.name())
            })?,
            None,
        );

        let sequence = Sequence::from(Vec::from(record.sequence()));

        writer
            .write_record(&fasta::Record::new(definition, sequence))
            .with_context(|| format!("Error writing FASTA record: {:?}", record))?;
    }

    Ok(())
}
