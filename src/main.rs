use clap::Parser;
use crossbeam::channel::{self, Receiver, Sender};
use flate2::{read::GzDecoder, write::GzEncoder, Compression};
use std::{
    collections::HashMap,
    fs::{self},
    io::{self, prelude::*, BufReader},
    sync::Arc,
    thread,
    time::{Duration, Instant},
};

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

#[derive(Parser, Debug)]
#[command(
    version,
    about = "Extract reads from a FASTQ file based on taxonomic classification via Kraken2."
)]
struct Args {
    #[arg(short, long)]
    kraken: String,
    #[arg(short, long)]
    taxid: i32,
    #[arg(short, long)]
    report: Option<String>,
    #[arg(long)]
    fastq: String,
    #[arg(long)]
    fastq2: Option<String>,
    #[arg(long)]
    output: String,
    #[arg(long)]
    output2: Option<String>,
    #[arg(long, default_value = "fast")]
    compression_mode: Option<String>,
    #[arg(long, action)]
    parents: bool,
    #[arg(long, action)]
    children: bool,
    #[arg(long)]
    no_compress: bool,
    #[arg(long)]
    exclude: bool,
    #[arg(long)]
    output_fasta: bool,
}

struct OutputConfig {
    no_compress: bool,
    compression_mode: Compression,
    output1: String,
    output2: Option<String>,
    output_fasta: bool,
}

impl OutputConfig {
    fn new(
        no_compress: bool,
        compression_mode: Option<String>,
        output1: String,
        output2: Option<String>,
        output_fasta: bool,
    ) -> Self {
        //detect compression mode
        let compression_mode = match compression_mode.as_deref() {
            // change compression mode to level(1-9) if specified
            Some("fast") => Compression::fast(),
            Some("default") => Compression::default(),
            Some("best") => Compression::best(),
            _ => {
                eprintln!("Invalid compression mode. Using default compression.");
                Compression::default()
            }
        };
        OutputConfig {
            no_compress,
            compression_mode,
            output1,
            output2,
            output_fasta,
        }
    }
}

/// Reads a FASTQ file from the specified path and returns a buffered reader.
///
/// This function reads a FASTQ file from the given path and returns a buffered reader
/// It automatically handles both plain text and gzipped files based on the file's magic number.
///
/// # Arguments
///
/// * `path` - A string containing the path to the FASTQ file.
///
/// # Returns
///
/// A buffered reader containing the contents of the FASTQ file. The reader may be a
/// plain text reader or a gzip decompressor, depending on the file format.
fn read_fastq(path: &str) -> BufReader<Box<dyn io::Read + Send>> {
    // The gzip magic number is 0x1f8b
    const GZIP_MAGIC_NUMBER: [u8; 2] = [0x1f, 0x8b];

    let mut file = if let Ok(file) = fs::File::open(path) {
        file
    } else {
        panic!("Error opening or reading the file");
    };

    // Buffer to store the first two bytes of the file
    let mut buffer = [0u8; 2];

    if let Ok(size) = file.read(&mut buffer) {
        // reset the pointer back to the start
        file.seek(io::SeekFrom::Start(0)).ok();
        // check if the first two bytes match the gzip magic number => file is gzipped
        // otherwise, it's a plain text file
        if size == 2 && buffer == GZIP_MAGIC_NUMBER {
            let gzip_reader: Box<dyn io::Read + Send> = Box::new(GzDecoder::new(file));
            io::BufReader::new(gzip_reader)
        } else {
            let plain_reader: Box<dyn io::Read + Send> = Box::new(file);
            io::BufReader::new(plain_reader)
        }
    } else {
        panic!("Error reading from the file");
    }
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
fn process_kraken_output_line(kraken_output: &str) -> (i32, String) {
    let fields: Vec<&str> = kraken_output.split('\t').collect();
    let taxon_id = fields[2].parse::<i32>().expect("Error parsing taxon ID");
    let read_id = fields[1].to_string();
    (taxon_id, read_id)
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
) -> Arc<HashMap<String, i32>> {
    let mut reads_to_save = HashMap::new();
    let kraken_file = fs::File::open(kraken_path).expect("Error reading kraken output file");
    let reader = io::BufReader::new(kraken_file);
    let mut total_reads = 0;

    print!("  Processing kraken output...");
    io::stdout().flush().unwrap();
    for line_result in reader.lines() {
        let line = line_result.expect("Error reading kraken output line");
        let (taxon_id, read_id) = process_kraken_output_line(&line);
        if exclude {
            if !taxon_ids_to_save.contains(&taxon_id) {
                reads_to_save.insert(read_id.clone(), taxon_id);
            }
        } else if taxon_ids_to_save.contains(&taxon_id) {
            reads_to_save.insert(read_id.clone(), taxon_id);
        }
        total_reads += 1;
    }
    println!("Done!");
    println!("  {} taxon IDs identified", taxon_ids_to_save.len());
    println!(
        "  {} total reads | {} reads to save.",
        total_reads,
        reads_to_save.len()
    );
    let reads_to_save: Arc<HashMap<String, i32>> = Arc::new(reads_to_save);
    reads_to_save
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
fn process_kraken_report_line(kraken_report: &str) -> (i32, i32) {
    let fields: Vec<&str> = kraken_report.split('\t').collect();
    let taxon_id = fields[4]
        .parse::<i32>()
        .expect("Error parsing taxon ID in kraken report");
    let mut spaces = 0;

    for char in fields[5].chars() {
        if char == ' ' {
            spaces += 1;
        } else {
            break;
        }
    }

    let level = spaces / 2;
    (taxon_id, level)
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
) -> (Vec<Tree>, HashMap<i32, usize>) {
    // will store the tree
    let mut nodes = Vec::new();
    // taxonid -> index in the nodes vector
    let mut taxon_map = HashMap::new();

    let report_file = if let Ok(report_file) = fs::File::open(report_path) {
        report_file
    } else {
        panic!("Error opening or reading the file");
    };
    {
        let reader = io::BufReader::new(report_file);
        let mut prev_index = None;

        for line in reader.lines().flatten() {
            let (taxon_id, level_num) = process_kraken_report_line(&line);
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
    (nodes, taxon_map)
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
fn collect_taxons_to_save(args: &Args) -> Vec<i32> {
    let mut taxon_ids_to_save = Vec::new();
    if args.report.is_some() {
        print!("  Processing kraken report...");
        io::stdout().flush().unwrap();
        let report_path = args.report.clone().unwrap();
        let (nodes, taxon_map) = build_tree_from_kraken_report(args.taxid, report_path);
        if args.children {
            let mut children = Vec::new();
            extract_children(&nodes, taxon_map[&args.taxid], &mut children);
            taxon_ids_to_save.extend(&children);
        } else if args.parents {
            taxon_ids_to_save.extend(extract_parents(&taxon_map, &nodes, args.taxid));
        } else {
            taxon_ids_to_save.push(args.taxid);
        }
        println!("Done!");
    } else {
        taxon_ids_to_save.push(args.taxid);
    }
    taxon_ids_to_save
}

/// Parse a FASTQ file and send reads to writer thread.
///
/// This function reads a fastq file and extracts read IDs and sequences.
/// It compares the read IDs against a given HashMap of read IDs (`reads_to_save`) and sends
/// the sequences of matching read IDs to the writer thread.
///
/// # Arguments
///
/// * `in_buf` - A buffered reader for the Fastq file.
/// * `tx` - A channel sender for sending sequences of selected reads.
/// * `reads_to_save` - A HashMap containing read IDs and corresponding taxon IDs.
fn parse_fastq(
    in_buf: io::BufReader<Box<dyn Read + Send>>,
    tx: &Sender<Vec<u8>>,
    reads_to_save: Arc<HashMap<String, i32>>,
    output_fasta: bool,
) {
    let mut num_lines = 0;
    let mut num_reads = 0;
    let mut current_id: String = String::new();
    //let mut stdout = BufWriter::new(io::stdout().lock());
    let mut line_bytes = Vec::new();
    let mut last_progress_update = Instant::now(); // For throttling progress updates
    const PROGRESS_UPDATE_INTERVAL: Duration = Duration::from_millis(1000);

    println!("  Reading fastq:");

    for line in in_buf.lines() {
        let line = line.expect("Error reading fastq line");
        line_bytes.clear();
        line_bytes.extend_from_slice(line.as_bytes());
        num_lines += 1;

        if num_lines % 4 == 1 {
            let fields: Vec<&str> = line.split_whitespace().collect();
            let read_id = &fields[0][1..];
            current_id = read_id.to_string();
            num_reads += 1;
        }

        if reads_to_save.contains_key(&current_id) {
            // if we want to outout a fasta, and are on first line of read - the ID
            if output_fasta && num_lines % 4 == 1 {
                let fasta_header = format!("> {}", current_id);
                tx.send(fasta_header.as_bytes().to_vec()).unwrap();
            // if we want to output a fasta, and are on the second line of the read - the sequence
            } else if output_fasta && num_lines % 4 == 2 {
                tx.send(line_bytes.to_vec()).unwrap();
            } else if !output_fasta {
                tx.send(line_bytes.to_vec()).unwrap();
            }
        }
        // Throttle progress updates - need to fix for multi-threading
        if last_progress_update.elapsed() >= PROGRESS_UPDATE_INTERVAL {
            print!("\rProcessed {} reads", num_reads);
            last_progress_update = Instant::now();
            io::stdout().flush().unwrap();
        }
    }
}

fn write_output_file(
    out_file: fs::File,
    rx: Receiver<Vec<u8>>,
    compression_mode: Compression,
    no_compress: bool,
) {
    let mut out_buf: Box<dyn io::Write> = if no_compress {
        Box::new(io::BufWriter::new(out_file))
    } else {
        Box::new(io::BufWriter::new(GzEncoder::new(
            out_file,
            compression_mode,
        )))
    };

    for data in rx {
        out_buf
            .write_all(&data)
            .and_then(|_| out_buf.write_all(b"\n"))
            .expect("Error writing to output file");
    }
    out_buf.flush().expect("Error flushing output buffer");
}

fn process_single_end(
    output_config: OutputConfig,
    reads_to_save_arc: Arc<HashMap<String, i32>>,
    in_buf: io::BufReader<Box<dyn Read + Send>>,
) {
    let (tx, rx) = channel::unbounded::<Vec<u8>>();
    let writer_thread = thread::spawn(move || {
        let out_file = fs::File::create(output_config.output1).expect("Error creating output file");
        write_output_file(
            out_file,
            rx,
            output_config.compression_mode,
            output_config.no_compress,
        );
    });
    let reader_thread = thread::spawn({
        let reads_to_save_arc = reads_to_save_arc.clone();
        move || {
            parse_fastq(in_buf, &tx, reads_to_save_arc, output_config.output_fasta);
        }
    });
    println!("  Processing is done. Writing is in progress...");
    writer_thread.join().unwrap();
    reader_thread.join().unwrap();
}

fn process_paired_end(
    output_config: OutputConfig,
    reads_to_save_arc: Arc<HashMap<String, i32>>,
    in_buf1: io::BufReader<Box<dyn Read + Send>>,
    in_buf2: io::BufReader<Box<dyn Read + Send>>,
) {
    let (tx1, rx1) = channel::unbounded::<Vec<u8>>();
    let (tx2, rx2) = channel::unbounded::<Vec<u8>>();

    let writer_thread1 = thread::spawn(move || {
        let out_file = fs::File::create(output_config.output1).expect("Error creating output file");
        write_output_file(
            out_file,
            rx1,
            output_config.compression_mode,
            output_config.no_compress,
        );
    });

    let writer_thread2 = thread::spawn(move || {
        let out_file =
            fs::File::create(output_config.output2.unwrap()).expect("Error creating output file");
        write_output_file(
            out_file,
            rx2,
            output_config.compression_mode,
            output_config.no_compress,
        );
    });
    let reader_thread1 = thread::spawn({
        let reads_to_save_arc = reads_to_save_arc.clone();
        move || {
            parse_fastq(in_buf1, &tx1, reads_to_save_arc, output_config.output_fasta);
        }
    });
    let reader_thread2 = thread::spawn({
        let reads_to_save_arc = reads_to_save_arc.clone();
        move || {
            parse_fastq(in_buf2, &tx2, reads_to_save_arc, output_config.output_fasta);
        }
    });

    writer_thread1.join().unwrap();
    writer_thread2.join().unwrap();
    reader_thread1.join().unwrap();
    reader_thread2.join().unwrap();
}

fn main() {
    let args = Args::parse();
    //check if paired-end reads are provided
    let paired = args.fastq2.is_some();
    if paired && args.output2.is_none() {
        panic!("Error: Paired-end reads provided but no output2 file specified");
    }
    if args.output2.is_some() && !paired {
        panic!("Error: output2 file specified but no paired-end reads provided");
    }

    println!(">> Step 1: Collecting taxon and read IDs to save");
    let taxon_ids_to_save = collect_taxons_to_save(&args);
    io::stdout().flush().unwrap();
    let reads_to_save = process_kraken_output(args.kraken, args.exclude, taxon_ids_to_save);
    //let reads_to_save_arc: Arc<HashMap<String, i32>> = Arc::new(reads_to_save);

    //create output from struct
    let output_config = OutputConfig::new(
        args.no_compress,
        args.compression_mode,
        args.output,
        args.output2,
        args.output_fasta,
    );

    println!(">> Step 2: Extracting reads to save and creating output");
    if !paired {
        // we are single end
        let in_buf = read_fastq(&args.fastq);
        process_single_end(output_config, reads_to_save.clone(), in_buf);
    } else {
        //we are paired end and so we need to spawn two writer threads and two reader threads
        let in_buf1 = read_fastq(&args.fastq);
        let in_buf2 = read_fastq(&args.fastq2.unwrap());
        process_paired_end(output_config, reads_to_save.clone(), in_buf1, in_buf2);
    }

    println!("Writing complete.");
}
