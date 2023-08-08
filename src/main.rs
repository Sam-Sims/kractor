use clap::Parser;
use flate2::{read::GzDecoder, write::GzEncoder, Compression};
use std::{
    collections::HashMap,
    fs,
    io::{self, prelude::*},
    sync::mpsc::channel,
    thread,
};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[arg(short, long)]
    kraken: String,
    #[arg(short, long)]
    taxid: i32,
    // #[arg(short, long)]
    // report: String,
    #[arg(short, long)]
    fastq: String,
    #[arg(short, long)]
    output: String,
}

// fn process_kraken_report(kraken_report: &str) -> (i32, &str, i32) {
//     let fields: Vec<&str> = kraken_report.split('\t').collect();
//     let taxon_id = fields[4].parse::<i32>().unwrap();
//     let rank_code = fields[3];
//     let mut spaces = 0;

//     for char in fields[5].chars() {
//         if char == ' ' {
//             spaces = spaces + 1;
//         }
//         else {
//             break;
//         }
//     }

//     let level = spaces / 2;
//     return (taxon_id, rank_code, level);

// }

fn process_kraken_output(kraken_output: &str) -> (i32, String) {
    let fields: Vec<&str> = kraken_output.split('\t').collect();
    let taxon_id = fields[2].parse::<i32>().unwrap();
    let read_id = fields[1].to_string();
    return (taxon_id, read_id);
}

fn read_fastq(path: &str) -> io::BufReader<Box<dyn io::Read>> {
    const GZIP_MAGIC_NUMBER: [u8; 2] = [0x1f, 0x8b];

    let mut file = match fs::File::open(path) {
        Ok(file) => file,
        Err(_) => panic!("Error opening or reading the file"),
    };

    let mut buffer = [0u8; 2];
    if let Ok(size) = file.read(&mut buffer) {
        file.seek(io::SeekFrom::Start(0)).ok();
        if size == 2 && buffer == GZIP_MAGIC_NUMBER {
            println!("File is gzipped");
            let gzip_reader: Box<dyn io::Read> = Box::new(GzDecoder::new(file));
            io::BufReader::new(gzip_reader)
        } else {
            println!("File is not gzipped");
            let plain_reader: Box<dyn io::Read> = Box::new(file);
            io::BufReader::new(plain_reader)
        }
    } else {
        panic!("Error reading from the file");
    }
}

fn main() {
    let args = Args::parse();

    // let kraken_report = fs::read_to_string(args.report)
    //     .expect("Something went wrong reading the file");

    // for line in kraken_report.lines() {
    //     let line_values = process_kraken_report(line);
    //     let (taxon_id, rank_code,mut level) = line_values;
    // }
    let mut reads_to_save = HashMap::new();
    println!("Reading kraken output...");
    let kraken_output =
        fs::read_to_string(args.kraken).expect("Something went wrong reading the file");

    let total_reads = kraken_output.lines().count();
    for line in kraken_output.lines() {
        let line_values = process_kraken_output(line);
        let (taxon_id, read_id) = line_values;
        if taxon_id == args.taxid {
            reads_to_save.insert(read_id, taxon_id);
        }
    }
    println!(
        "Done! {} total reads | {} reads to save.",
        reads_to_save.len(),
        total_reads
    );

    let mut num_lines = 0;
    let mut num_reads = 0;
    let mut current_id: String = String::new();

    let in_buf = read_fastq(&args.fastq);

    let (tx, rx) = channel::<Vec<u8>>();
    let writer_thread = thread::spawn(move || {
        let out_file = fs::File::create(args.output).unwrap();
        let out_gzip = GzEncoder::new(out_file, Compression::fast());
        let mut out_buf = io::BufWriter::new(out_gzip);

        for data in rx {
            out_buf.write_all(&data);
            out_buf.write_all(b"\n");
            out_buf.flush();
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
            let data_to_write = line_bytes.to_vec();
            tx.send(data_to_write).unwrap();
            print!("Processed {} reads\r", num_reads);
            io::stdout().flush().unwrap();
        }
    }
    println!("Processing is done. Writing is in progress...");
    drop(tx);

    writer_thread.join().unwrap();

    println!("Writing complete.");
}
