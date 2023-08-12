use clap::Parser;
use flate2::{read::GzDecoder, write::GzEncoder, Compression};
use std::{
    collections::HashMap,
    fs::{self},
    io::{self, prelude::*},
    sync::mpsc::channel,
    thread,
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
    #[arg(short, long)]
    fastq: String,
    #[arg(short, long)]
    output: String,
    #[arg(long, default_value = "fast")]
    compression: Option<String>,
    #[arg(long, action)]
    parents: bool,
    #[arg(long, action)]
    children: bool,
}

fn process_kraken_output(kraken_output: &str) -> (i32, String) {
    let fields: Vec<&str> = kraken_output.split('\t').collect();
    let taxon_id = fields[2].parse::<i32>().expect("Error parsing taxon ID");
    let read_id = fields[1].to_string();
    return (taxon_id, read_id);
}

fn read_fastq(path: &str) -> io::BufReader<Box<dyn io::Read>> {
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
            let gzip_reader: Box<dyn io::Read> = Box::new(GzDecoder::new(file));
            io::BufReader::new(gzip_reader)
        } else {
            let plain_reader: Box<dyn io::Read> = Box::new(file);
            io::BufReader::new(plain_reader)
        }
    } else {
        panic!("Error reading from the file");
    }
}

fn extract_parents(
    taxon_to_index: &HashMap<i32, usize>,
    nodes: &Vec<Tree>,
    taxon_id: i32,
) -> Vec<i32> {
    // Backtracking traversal from the given taxon_id to the root
    let mut parents = Vec::new();
    parents.push(taxon_id);
    let mut curr_index = taxon_to_index[&taxon_id];

    while let Some(parent_index) = nodes[curr_index].parent {
        parents.push(nodes[parent_index].taxon_id);
        curr_index = parent_index;
    }

    parents
}

fn extract_children(nodes: &Vec<Tree>, start_index: usize, result: &mut Vec<i32>) {
    // recursive post-order traversal of the tree
    for &child_index in &nodes[start_index].children {
        extract_children(nodes, child_index, result);
    }
    result.push(nodes[start_index].taxon_id);
}

fn process_kraken_report_line(kraken_report: &str) -> (i32, i32) {
    let fields: Vec<&str> = kraken_report.split('\t').collect();
    let taxon_id = fields[4].parse::<i32>().unwrap();
    let mut spaces = 0;

    for char in fields[5].chars() {
        if char == ' ' {
            spaces = spaces + 1;
        } else {
            break;
        }
    }

    let level = spaces / 2;
    (taxon_id, level)
}

fn process_kraken_report(
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

        for line in reader.lines() {
            if let Ok(line) = line {
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

                // add the current node to the parent's children
                if let Some(parent_index) = prev_index {
                    nodes[parent_index].children.push(curr_index);
                }

                prev_index = Some(curr_index);

                // if the current taxon is the one we want to save, add it to the map
                if taxon_id == taxon_to_save {
                    taxon_map.insert(taxon_id, curr_index);
                }
            }
        }
    }
    (nodes, taxon_map)
}

fn main() {
    let args = Args::parse();

    let compression_mode = match args.compression.as_deref() {
        Some("fast") => Compression::fast(),
        Some("default") => Compression::default(),
        Some("best") => Compression::best(),
        _ => {
            eprintln!("Invalid compression mode. Using default compression.");
            Compression::default()
        }
    };
    // if --report is provided, read the report and build the tree
    let mut taxon_ids_to_save = Vec::new();

    if args.report.is_some() {
        println!("Reading kraken report");
        let (nodes, taxon_to_index) =
            process_kraken_report(args.taxid, args.report.unwrap().to_string());
        if args.children {
            // declearing children here so that it's not reallocated every time the function is called
            let mut children = Vec::new();
            extract_children(&nodes, taxon_to_index[&args.taxid], &mut children);
            //println!("{:?}", children);
            taxon_ids_to_save.extend(children);
        }
        if args.parents {
            //println!("{:?}", extract_parents(&taxon_to_index, &nodes, args.taxid));
            taxon_ids_to_save.extend(extract_parents(&taxon_to_index, &nodes, args.taxid));
        }
    } else {
        taxon_ids_to_save.push(args.taxid);
    }

    println!("Saving reads from taxon(s): {:?}", taxon_ids_to_save);

    let mut reads_to_save = HashMap::new();
    println!("Reading kraken output");
    let kraken_output = fs::read_to_string(args.kraken).expect("Error reading kraken output file");

    let total_reads = kraken_output.lines().count();
    for line in kraken_output.lines() {
        let (taxon_id, read_id) = process_kraken_output(line);
        if taxon_ids_to_save.contains(&taxon_id) {
            reads_to_save.insert(read_id, taxon_id);
        }
    }
    println!(
        "Done! {} total reads | {} reads to save.",
        total_reads,
        reads_to_save.len()
    );

    let mut num_lines = 0;
    let mut num_reads = 0;
    let mut current_id: String = String::new();

    let in_buf = read_fastq(&args.fastq);

    let (tx, rx) = channel::<Vec<u8>>();
    let writer_thread = thread::spawn(move || {
        let out_file = fs::File::create(args.output).expect("Error creating output file");
        let out_gzip = GzEncoder::new(out_file, compression_mode);
        let mut out_buf = io::BufWriter::new(out_gzip);

        for data in rx {
            out_buf
                .write_all(&data)
                .and_then(|_| out_buf.write_all(b"\n"))
                .and_then(|_| out_buf.flush())
                .expect("Error writing to output file");
        }
    });
    println!("Reading fastq:");

    for line in in_buf.lines() {
        let line = line.unwrap();
        let line_bytes = line.as_bytes();
        num_lines += 1;

        match num_lines % 4 {
            1 => {
                let fields: Vec<&str> = line.split_whitespace().collect();
                let read_id = fields[0].to_string();
                current_id = read_id[1..].to_string();
                num_reads += 1;
            }
            _ => {}
        };

        if reads_to_save.contains_key(&current_id) {
            //let data_to_write = line_bytes.to_vec();
            tx.send(line_bytes.to_vec()).unwrap();
            print!("Processed {} reads\r", num_reads);
            io::stdout().flush().unwrap();
        }
    }
    println!("Processing is done. Writing is in progress...");
    drop(tx);

    writer_thread.join().unwrap();

    println!("Writing complete.");
}
