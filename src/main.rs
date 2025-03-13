pub use crate::cli::Cli;
use chrono::Local;
use clap::Parser;
use crossbeam::channel::{self, Receiver, Sender};
use env_logger::{fmt::Color, Builder};
use log::{debug, info, trace, LevelFilter};
use noodles::{
    fasta::{
        self,
        record::{Definition, Sequence},
    },
    fastq,
};
use std::{
    collections::{HashMap, HashSet},
    fs::{self},
    io::{self, prelude::*, BufReader},
    path::Path,
    sync::{Arc, Mutex},
    thread,
    time::{Duration, Instant},
};
use lazy_static::lazy_static;
use anyhow::{Result, Context, anyhow, bail};

mod cli;

#[derive(Debug, Clone)]
struct Tree {
    taxon_id: i32,
    level_num: i32,
    children: Vec<usize>,
    parent: Option<usize>,
}

impl Tree {
    fn new(taxon_id: i32, level_num: i32, parent: Option<usize>) -> Tree {
        Tree {
            taxon_id,
            level_num,
            children: Vec::new(),
            parent,
        }
    }
}

lazy_static! {
    static ref TAXON_ID_COUNT: Arc<Mutex<usize>> = Arc::new(Mutex::new(0));
    static ref TAXON_IDS: Arc<Mutex<Vec<i32>>> = Arc::new(Mutex::new(Vec::new()));
    static ref TOTAL_READS: Arc<Mutex<usize>> = Arc::new(Mutex::new(0));
    static ref READS_TO_EXTRACT: Arc<Mutex<usize>> = Arc::new(Mutex::new(0));
}

/// Parses a Kraken output line to extract taxon ID and read ID.
///
/// This function takes a Kraken output line and processes it to extract the taxon ID and read ID.
/// It is called in the process_kraken_output function.
///
/// # Arguments
///
/// * `kraken_output` - A string representing a single line from Kraken output.
///
/// # Returns
///
/// A tuple containing the extracted taxon ID and the read ID
fn process_kraken_output_line(kraken_output: &str) -> Result<(i32, String)> {
    let fields: Vec<&str> = kraken_output.split('\t').collect();
    if fields.len() < 4 {
        bail!("Invalid kraken output format: Expected at least 4 tab-separated fields, got {}. Line: '{}'",
              fields.len(), kraken_output);
    }
    let taxon_id = fields[2].parse::<i32>().context("Error parsing taxon ID")?;
    let read_id = fields[1].to_string();
    Ok((taxon_id, read_id))
}

/// Processes the Kraken output file to extract read ID
///
/// This function takes the kraken output file and processes each line to extract the taxon id and read id.
/// Read IDs that match the taxon IDs to save are stored in a hashmap.
///
/// # Arguments
///
/// `kraken_path` - A string containing the path to the Kraken output file.
/// `exclude` - A boolean indicating whether to exclude or include the taxon IDs to save.
/// `taxon_ids_to_save` - A vector containing the taxon IDs to save.
///
/// # Returns
///
/// A hashmap containing the read IDs to save as keys and the taxon IDs as values.
fn process_kraken_output(
    kraken_path: String,
    exclude: bool,
    taxon_ids_to_save: Vec<i32>,
) -> Result<Arc<HashSet<String>>> {
    info!("Processing kraken output...");
    let mut reads_to_save = HashSet::new();
    let kraken_file = fs::File::open(&kraken_path)
        .with_context(|| format!("Failed to open kraken output file: {}", kraken_path))?;
    let reader = io::BufReader::new(kraken_file);
    let mut total_reads = 0;

    for line_result in reader.lines() {
        let line = line_result.context("Error reading kraken output line")?;
        let (taxon_id, read_id) = process_kraken_output_line(&line)?;
        if exclude {
            if !taxon_ids_to_save.contains(&taxon_id) {
                reads_to_save.insert(read_id);
            }
        } else if taxon_ids_to_save.contains(&taxon_id) {
            reads_to_save.insert(read_id);
        }
        total_reads += 1;
    }
    let mut taxon_id_count = TAXON_ID_COUNT.lock()
        .map_err(|e| anyhow!("Failed to lock TAXON_ID_COUNT mutex: {}", e))?;
    let mut taxon_ids = TAXON_IDS.lock()
        .map_err(|e| anyhow!("Failed to lock TAXON_IDS mutex: {}", e))?;
    let mut total_read_count = TOTAL_READS.lock()
        .map_err(|e| anyhow!("Failed to lock TOTAL_READS mutex: {}", e))?;
    let mut reads_to_extract = READS_TO_EXTRACT.lock()
        .map_err(|e| anyhow!("Failed to lock READS_TO_EXTRACT mutex: {}", e))?;

    *taxon_id_count = taxon_ids_to_save.len();
    *taxon_ids = taxon_ids_to_save;
    *total_read_count = total_reads; // Update this with the actual total_reads value
    *reads_to_extract = reads_to_save.len();

    let reads_to_save: Arc<HashSet<String>> = Arc::new(reads_to_save);
    Ok(reads_to_save)
}

/// Extracts the taxon ID of all parents for a given taxon ID.
///
/// This function implements a backtracking traversal from the specified `taxon_id` to the root.
///
/// # Arguments
///
/// * `taxon_map` - Mapping of taxon IDs to their corresponding indices in the `nodes` vector.
/// * `nodes` - The tree.
/// * `taxon_id` - The taxon ID for which to extract the lineage of parent taxon IDs.
///
/// # Returns
///
/// A vector containing the taxon IDs of the lineage of parent nodes, including the provided taxon ID.
fn extract_parents(taxon_map: &HashMap<i32, usize>, nodes: &[Tree], taxon_id: i32) -> Vec<i32> {
    // Backtracking traversal from the given taxon_id to the root
    let mut parents = Vec::new();
    parents.push(taxon_id);
    let mut curr_index = taxon_map[&taxon_id];

    while let Some(parent_index) = nodes[curr_index].parent {
        parents.push(nodes[parent_index].taxon_id);
        curr_index = parent_index;
    }

    parents
}

/// Extracts the taxon IDs of children nodes from a given taxon ID.
///
/// This function implements a recursive post-order traversal of the tree starting from
/// the specified taxon. It collects the taxon IDs of child nodes and appends them to the
/// provided result vector.
///
/// # Arguments
///
/// * `nodes` - The tree.
/// * `start_index` - The node to start the traversal from.
/// * `result` - Stores the extracted child taxon IDs.
///
/// # Returns
///
/// A vector containing the taxon IDs of the children of the specified taxon ID, including the provided taxon ID.
fn extract_children(nodes: &Vec<Tree>, start_index: usize, result: &mut Vec<i32>) {
    // recursive post-order traversal of the tree
    for &child_index in &nodes[start_index].children {
        extract_children(nodes, child_index, result);
    }
    result.push(nodes[start_index].taxon_id);
}

/// Parses a Kraken report line to extract taxon ID and its corresponding level.
///
/// This function takes a Kraken report line and processes it to extract the taxon ID
/// and its taxonomic level. The level is calculated based on the indentation of the taxon name field.
/// It is called in the process_kraken_report function.
///
/// # Arguments
///
/// * `kraken_report` - A string representing a single line from a Kraken report.
///
/// # Returns
///
/// A tuple containing the extracted taxon ID and its corresponding level.
fn process_kraken_report_line(kraken_report: &str) -> Result<(i32, i32)> {
    let fields: Vec<&str> = kraken_report.split('\t').collect();
    if fields.len() < 6 {
        bail!("Invalid kraken report line format: Expected at least 6 tab-separated fields, got {}. Line: '{}'",
              fields.len(), kraken_report);
    }
    let taxon_id = fields[4]
        .parse::<i32>()
        .context("Error parsing taxon ID in kraken report")?;
    let mut spaces = 0;

    for char in fields[5].chars() {
        if char == ' ' {
            spaces += 1;
        } else {
            break;
        }
    }

    let level = spaces / 2;
    Ok((taxon_id, level))
}

/// Processes a Kraken report to build a tree of all taxa in the kraken report.
///
/// This function reads a Kraken report from the specified path and processes it to
/// construct a taxonomic tree. Each node corresponds to a taxon.
///
/// # Arguments
///
/// * `taxon_to_save` - The taxon ID that needs to be saved.
/// * `report_path` - A string containing the path to the Kraken report file.
///
/// # Returns
///
/// A tuple containing the tree and a hashmap mapping the saved taxon IDs to the tree.
fn build_tree_from_kraken_report(
    taxon_to_save: i32,
    report_path: String,
) -> Result<(Vec<Tree>, HashMap<i32, usize>)> {
    debug!("Building taxonomic tree from kraken report");
    // will store the tree
    let mut nodes = Vec::new();
    // taxonid -> index in the nodes vector
    let mut taxon_map = HashMap::new();

    let report_file = fs::File::open(&report_path)
        .with_context(|| format!("Failed to open kraken report file: {}", report_path))?;

    {
        let reader = BufReader::new(report_file);
        let mut prev_index = None;

        for line in reader.lines().flatten() {
            let (taxon_id, level_num) = process_kraken_report_line(&line)?;
            // if taxon_id == 0, it's an unclassified read so we can skip
            if taxon_id == 0 {
                continue;
            }
            // 1 will be the root of the tree
            if taxon_id == 1 {
                let root_node = Tree::new(taxon_id, level_num, None);
                prev_index = Some(nodes.len());
                nodes.push(root_node);
            }
            // if the current level is not the same as the previous level + 1, then we are not at the correct parent, and need to move up the tree
            while let Some(parent_index) = prev_index {
                if level_num != nodes[parent_index].level_num + 1 {
                    prev_index = nodes[parent_index].parent;
                } else {
                    break;
                }
            }
            // once we have the correct parent, we can add the current node to the tree
            let curr_node = Tree::new(taxon_id, level_num, prev_index);
            let curr_index = nodes.len();
            nodes.push(curr_node);

            // add the current node
            if let Some(parent_index) = prev_index {
                nodes[parent_index].children.push(curr_index);
            }

            prev_index = Some(curr_index);

            // if the current taxon is one we want to save, add it to the map
            if taxon_id == taxon_to_save {
                taxon_map.insert(taxon_id, curr_index);
            }
        }
    }

    if !taxon_map.contains_key(&taxon_to_save) {
        bail!("Taxon ID {} not found in the kraken report", taxon_to_save);
    }

    Ok((nodes, taxon_map))
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
fn collect_taxons_to_save(args: &Cli) -> Result<Vec<i32>> {
    //TODO refactor this not to use the entire args struct
    let mut taxon_ids_to_save = Vec::new();
    if args.report.is_some() {
        info!("Processing kraken report...");
        let report_path = args.report.clone().unwrap();
        let (nodes, taxon_map) = build_tree_from_kraken_report(args.taxid, report_path)?;
        if args.children {
            debug!("Extracting children");
            let mut children = Vec::new();
            extract_children(&nodes, taxon_map[&args.taxid], &mut children);
            taxon_ids_to_save.extend(&children);
        } else if args.parents {
            debug!("Extracting parents");
            taxon_ids_to_save.extend(extract_parents(&taxon_map, &nodes, args.taxid));
        } else {
            taxon_ids_to_save.push(args.taxid);
        }
    } else {
        debug!(
            "No kraken report provided - extracting reads for taxon ID {} only",
            args.taxid
        );
        taxon_ids_to_save.push(args.taxid);
    }
    debug!("Taxon IDs identified: {:?}", taxon_ids_to_save);
    if taxon_ids_to_save.is_empty() {
        bail!("No taxon IDs were identified for extraction");
    }
    Ok(taxon_ids_to_save)
}


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
fn parse_fastq(file_path: &str, reads_to_save: Arc<HashSet<String>>, tx: &Sender<fastq::Record>) -> Result<()> {
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
        let read_name_string = std::str::from_utf8(read_name_bytes)
            .context("Invalid UTF-8 sequence in read name")?;

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
fn infer_compression(file_path: &String) -> niffler::compression::Format {
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
fn write_output_fastq(
    rx: Receiver<fastq::Record>,
    out_file: String,
    output_type: Option<niffler::Format>,
    compression_level: niffler::Level,
) -> Result<()> {
    let compression_type = match output_type {
        Some(output_type) => {
            debug!("Output type overridden as: {:?}", output_type);
            output_type
        }
        None => {
            let inferred_type = infer_compression(&out_file);
            debug!("Inferred output compression type as: {:?}", inferred_type);
            inferred_type
        }
    };

    debug!(
        "Output compression level specified as: {:?}",
        compression_level
    );
    debug!("Creating output file: {:?}", out_file);

    let out_file_path = out_file.clone();
    let out_file = fs::File::create(&out_file)
        .with_context(|| format!("Failed to create output file: {}", out_file_path))?;

    let file_handle = Box::new(io::BufWriter::new(out_file));
    let writer = niffler::get_writer(file_handle, compression_type, compression_level)
        .context("Failed to create compressed writer")?;

    let mut fastq_writer = fastq::Writer::new(writer);

    for record in rx {
        fastq_writer.write_record(&record)
            .context("Error writing FASTQ record")?;
    }

    Ok(())
}


fn write_output_fasta(rx: Receiver<fastq::Record>, out_file: String) -> Result<()> {
    debug!("Creating output file: {:?}", out_file);

    let out_file_path = out_file.clone();
    let out_file = fs::File::create(&out_file)
        .with_context(|| format!("Failed to create output file: {}", out_file_path))?;

    let mut writer = fasta::Writer::new(out_file);

    for record in rx {
        let definition = Definition::new(
            std::str::from_utf8(record.name())
                .context("Invalid UTF-8 sequence in read name")?,
            None
        );

        let sequence = Sequence::from(Vec::from(record.sequence()));

        writer.write_record(&fasta::Record::new(definition, sequence))
            .context("Error writing FASTA record")?;
    }

    Ok(())
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
    reads_to_save: Arc<HashSet<String>>,
    input: Vec<String>,
    output: Vec<String>,
    compression_type: Option<niffler::Format>,
    compression_level: niffler::Level,
    fasta: bool,
) -> Result<()> {
    let (tx, rx) = channel::unbounded::<fastq::Record>();
    let output_file = output[0].clone();

    let reader_thread = thread::spawn({
        trace!("Spawning reader thread");
        let reads_to_save_arc = reads_to_save.clone();
        let tx_clone = tx.clone();
        move || -> Result<()> {
            parse_fastq(&input[0], reads_to_save_arc, &tx_clone)
        }
    });

    let writer_thread = thread::spawn({
        trace!("Spawning writer thread");
        move || -> Result<()> {
            if !fasta {
                write_output_fastq(rx, output_file, compression_type, compression_level)
            } else {
                write_output_fasta(rx, output_file)
            }
        }
    });

    let reader_result = reader_thread.join()
        .map_err(|e| anyhow!("Reader thread panicked: {:?}", e))?;
    reader_result.context("Reader thread operation failed")?;

    info!("Processing complete. Writing is in progress...");

    let writer_result = writer_thread.join()
        .map_err(|e| anyhow!("Writer thread panicked: {:?}", e))?;
    writer_result.context("Writer thread operation failed")?;

    info!("Writing complete.");
    Ok(())
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
    reads_to_save: Arc<HashSet<String>>,
    input: Vec<String>,
    output: Vec<String>,
    compression_type: Option<niffler::Format>,
    compression_level: niffler::Level,
    fasta: bool,
) -> Result<()> {
    let (tx1, rx1) = channel::unbounded::<fastq::Record>();
    let (tx2, rx2) = channel::unbounded::<fastq::Record>();
    let input_file1 = input[0].clone();
    let input_file2 = input[1].clone();
    let output_file1 = output[0].clone();
    let output_file2 = output[1].clone();

    let reader_thread1 = thread::spawn({
        trace!("Spawning reader thread 1");
        let reads_to_save_arc = reads_to_save.clone();
        let tx_clone = tx1.clone();
        move || -> Result<()> {
            parse_fastq(&input_file1, reads_to_save_arc, &tx_clone)
        }
    });

    let reader_thread2 = thread::spawn({
        trace!("Spawning reader thread 2");
        let reads_to_save_arc = reads_to_save.clone();
        let tx_clone = tx2.clone();
        move || -> Result<()> {
            parse_fastq(&input_file2, reads_to_save_arc, &tx_clone)
        }
    });

    let writer_thread1 = thread::spawn({
        trace!("Spawning writer thread 1");
        move || -> Result<()> {
            if !fasta {
                write_output_fastq(rx1, output_file1, compression_type, compression_level)
            } else {
                write_output_fasta(rx1, output_file1)
            }
        }
    });

    let writer_thread2 = thread::spawn({
        trace!("Spawning writer thread 2");
        move || -> Result<()> {
            if !fasta {
                write_output_fastq(rx2, output_file2, compression_type, compression_level)
            } else {
                write_output_fasta(rx2, output_file2)
            }
        }
    });

    let reader_result1 = reader_thread1.join()
        .map_err(|e| anyhow!("Reader thread 1 panicked: {:?}", e))?;
    reader_result1.context("Reader thread 1 operation failed")?;

    let reader_result2 = reader_thread2.join()
        .map_err(|e| anyhow!("Reader thread 2 panicked: {:?}", e))?;
    reader_result2.context("Reader thread 2 operation failed")?;

    info!("Processing is done. Writing is in progress...");

    let writer_result1 = writer_thread1.join()
        .map_err(|e| anyhow!("Writer thread 1 panicked: {:?}", e))?;
    writer_result1.context("Writer thread 1 operation failed")?;

    let writer_result2 = writer_thread2.join()
        .map_err(|e| anyhow!("Writer thread 2 panicked: {:?}", e))?;
    writer_result2.context("Writer thread 2 operation failed")?;

    Ok(())
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
    let taxon_ids_to_save = collect_taxons_to_save(&args)?;
    let reads_to_save = process_kraken_output(args.kraken, args.exclude, taxon_ids_to_save)?;

    //check if paired-end reads are provided
    let paired = args.input.len() == 2;
    if paired {
        info!("Detected two input files. Assuming paired-end reads");
        process_paired_end(
            reads_to_save.clone(),
            args.input,
            args.output,
            args.output_type,
            args.compression_level,
            args.output_fasta,
        )?;
    } else {
        info!("Detected one input file. Assuming single-end reads");
        process_single_end(
            reads_to_save.clone(),
            args.input,
            args.output,
            args.output_type,
            args.compression_level,
            args.output_fasta,
        )?;
    }

    info!("Complete!");

    if !args.no_json {
        let stats = serde_json::json!({
            "taxon_id_count": *TAXON_ID_COUNT.lock()
                .map_err(|e| anyhow!("Failed to lock TAXON_ID_COUNT mutex: {}", e))?,
            "taxon_ids": *TAXON_IDS.lock()
                .map_err(|e| anyhow!("Failed to lock TAXON_IDS mutex: {}", e))?,
            "reads_in": *TOTAL_READS.lock()
                .map_err(|e| anyhow!("Failed to lock TOTAL_READS mutex: {}", e))?,
            "reads_out": *READS_TO_EXTRACT.lock()
                .map_err(|e| anyhow!("Failed to lock READS_TO_EXTRACT mutex: {}", e))?,
            "input_format": if paired { "paired-end" } else { "single-end" },
            "output_format": if args.output_fasta { "fasta" } else { "fastq" },
        });

        let stats_json = serde_json::to_string_pretty(&stats)
            .context("Failed to serialize stats to JSON")?;

        println!("{}", stats_json);
    }

    Ok(())
}
