use color_eyre::eyre::{Context, Result};
use crossbeam::channel::{Receiver, Sender};
use fxhash::FxHashSet;
use log::{debug, trace};
use noodles::fasta::record::{Definition, Sequence};
use noodles::{fasta, fastq};
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};
use std::{fs, io};

pub fn parse_fastq(
    file_path: &PathBuf,
    reads_to_save: &FxHashSet<Vec<u8>>,
    tx: &Sender<fastq::Record>,
) -> Result<usize> {
    const PROGRESS_UPDATE_INTERVAL: Duration = Duration::from_millis(1500);
    
    let mut num_reads = 0;
    let mut last_progress_update = Instant::now();

    let (reader, format) = niffler::from_path(file_path)
        .wrap_err_with(|| format!("Failed to open fastq file: {}", file_path.display()))?;
    debug!(
        "Detected input compression type for file {} as: {format:?}",
        file_path.display()
    );
    let reader = BufReader::new(reader);
    let mut fastq_reader = fastq::Reader::new(reader);

    for (record_idx, result) in fastq_reader.records().enumerate() {
        let record = result
            .wrap_err_with(|| format!("Error reading FASTQ record at position {record_idx}"))?;

        let read_id = record.name();
        if reads_to_save.contains(&read_id.to_vec()) {
            tx.send(record).wrap_err("Error sending record")?;
        }
        num_reads += 1;

        if last_progress_update.elapsed() >= PROGRESS_UPDATE_INTERVAL {
            trace!("Processed {num_reads} reads");
            last_progress_update = Instant::now();
        }
    }

    Ok(num_reads)
}

fn infer_compression(file_path: &PathBuf) -> niffler::compression::Format {
    let path = Path::new(&file_path);
    let ext = path.extension().unwrap().to_str().unwrap();
    match ext {
        "gz" => niffler::compression::Format::Gzip,
        "bz2" => niffler::compression::Format::Bzip,
        _ => niffler::compression::Format::No,
    }
}

pub fn write_output_fastq(
    rx: Receiver<fastq::Record>,
    out_file: &PathBuf,
    output_type: Option<niffler::Format>,
    compression_level: niffler::Level,
) -> Result<usize> {
    let mut read_output_count = 0;
    let compression_type = if let Some(output_type) = output_type {
        debug!("Output type overridden as: {output_type:?}");
        output_type
    } else {
        let inferred_type = infer_compression(out_file);
        debug!("Inferred output compression type as: {inferred_type:?}");
        inferred_type
    };

    debug!("Output compression level specified as: {compression_level:?}");
    debug!("Creating output file: {}", out_file.display());

    fs::create_dir_all(out_file.parent().unwrap())
        .wrap_err_with(|| format!("Failed to create output directory: {}", out_file.display()))?;

    let out_file = fs::File::create(out_file)
        .wrap_err_with(|| format!("Failed to create output file: {}", out_file.display()))?;

    let file_handle = Box::new(io::BufWriter::new(out_file));
    let writer = niffler::get_writer(file_handle, compression_type, compression_level)
        .wrap_err("Failed to create niffler writer")?;

    let mut fastq_writer = fastq::Writer::new(writer);

    for record in rx {
        fastq_writer
            .write_record(&record)
            .wrap_err_with(|| format!("Error writing FASTQ record: {record:?}"))?;
        read_output_count += 1;
    }

    Ok(read_output_count)
}

pub fn write_output_fasta(rx: Receiver<fastq::Record>, out_file: &PathBuf) -> Result<usize> {
    debug!("Creating output file: {}", out_file.display());
    let mut total_read_count = 0;
    let out_file = fs::File::create(out_file)
        .wrap_err_with(|| format!("Failed to create output file: {}", out_file.display()))?;

    let mut writer = fasta::Writer::new(out_file);

    for record in rx {
        let definition = Definition::new(
            std::str::from_utf8(record.name()).wrap_err_with(|| {
                format!("Invalid UTF-8 sequence in read name: {:?}", record.name())
            })?,
            None,
        );

        let sequence = Sequence::from(Vec::from(record.sequence()));

        writer
            .write_record(&fasta::Record::new(definition, sequence))
            .wrap_err_with(|| format!("Error writing FASTA record: {record:?}"))?;
        total_read_count += 1;
    }

    Ok(total_read_count)
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::fastq;
    use std::fs::File;
    use std::io::{Read, Write};
    use tempfile::tempdir;

    #[test]
    fn test_infer_compression_gzip() {
        let file_path = PathBuf::from("test.gz");
        let compression = infer_compression(&file_path);

        assert_eq!(compression, niffler::compression::Format::Gzip);
    }

    #[test]
    fn test_infer_compression_bzip() {
        let file_path = PathBuf::from("test.bz2");
        let compression = infer_compression(&file_path);

        assert_eq!(compression, niffler::compression::Format::Bzip);
    }

    #[test]
    fn test_infer_compression_no_compression() {
        let file_path = PathBuf::from("test.fastq");
        let compression = infer_compression(&file_path);

        assert_eq!(compression, niffler::compression::Format::No);
    }

    #[test]
    fn test_parse_fastq_with_matches() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("test.fastq");
        let test_data = b"@read1\nAAAA\n+\n!!!!\n@read2\nGGGG\n+\n!!!!\n@read3\nTTTT\n+\n!!!!\n";
        let mut file = File::create(&file_path).unwrap();
        file.write_all(test_data).unwrap();
        let mut reads_to_save = FxHashSet::default();
        reads_to_save.insert(b"read1".to_vec());
        reads_to_save.insert(b"read3".to_vec());
        let (tx, rx) = crossbeam::channel::unbounded();
        parse_fastq(&file_path, &reads_to_save, &tx).unwrap();
        drop(tx);
        let results: Vec<fastq::Record> = rx.iter().collect();

        assert_eq!(results.len(), 2);
        assert_eq!(results[0].name(), b"read1");
        assert_eq!(results[1].name(), b"read3");
        assert_eq!(results[0].sequence(), b"AAAA");
        assert_eq!(results[1].sequence(), b"TTTT");
    }

    #[test]
    fn test_parse_fastq_with_no_matches() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("test.fastq");
        let test_data = b"@read1\nAAAA\n+\n!!!!\n@read2\nGGGG\n+\n!!!!\n@read3\nTTTT\n+\n!!!!\n";
        let mut file = File::create(&file_path).unwrap();
        file.write_all(test_data).unwrap();
        let mut reads_to_save = FxHashSet::default();
        reads_to_save.insert(b"read4".to_vec());
        reads_to_save.insert(b"read5".to_vec());
        let (tx, rx) = crossbeam::channel::unbounded();
        parse_fastq(&file_path, &reads_to_save, &tx).unwrap();
        drop(tx);
        let results: Vec<fastq::Record> = rx.iter().collect();

        assert_eq!(results.len(), 0);
    }

    #[test]
    fn test_parse_fastq_file_not_found() {
        let file_path = PathBuf::from("idontexist.fastq");
        let reads_to_save = FxHashSet::default();
        let (tx, _rx) = crossbeam::channel::unbounded();
        let result = parse_fastq(&file_path, &reads_to_save, &tx);

        assert!(result.is_err());
    }

    #[test]
    fn test_write_output_fastq_non_compressed() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("output.fastq");
        let (tx, rx) = crossbeam::channel::unbounded();
        let record1 = fastq::Record::new(
            fastq::record::Definition::new("read1", "read1"),
            "AAAA",
            "!!!!",
        );
        let record2 = fastq::Record::new(
            fastq::record::Definition::new("read2", "read2"),
            "GGGG",
            "!!!!",
        );
        tx.send(record1).unwrap();
        tx.send(record2).unwrap();
        drop(tx);
        let read_count = write_output_fastq(
            rx,
            &file_path,
            Some(niffler::compression::Format::No),
            niffler::Level::One,
        )
        .unwrap();
        let file_content = fs::read_to_string(file_path).unwrap();

        assert_eq!(read_count, 2);
        assert!(file_content.contains("@read1"));
        assert!(file_content.contains("AAAA"));
        assert!(file_content.contains("@read2"));
        assert!(file_content.contains("GGGG"));
    }

    #[test]
    fn test_write_output_fastq_gzip() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("output.fastq");
        let (tx, rx) = crossbeam::channel::unbounded();
        let record1 = fastq::Record::new(
            fastq::record::Definition::new("read1", "read1"),
            "AAAA",
            "!!!!",
        );
        let record2 = fastq::Record::new(
            fastq::record::Definition::new("read2", "read2"),
            "GGGG",
            "!!!!",
        );
        tx.send(record1).unwrap();
        tx.send(record2).unwrap();
        drop(tx);
        let read_count = write_output_fastq(
            rx,
            &file_path,
            Some(niffler::compression::Format::Gzip),
            niffler::Level::One,
        )
        .unwrap();
        let reader =
            niffler::get_reader(Box::new(BufReader::new(File::open(&file_path).unwrap()))).unwrap();
        let mut decompressed = String::new();
        let mut decompressed_reader = BufReader::new(reader.0);
        decompressed_reader
            .read_to_string(&mut decompressed)
            .unwrap();

        assert_eq!(read_count, 2);
        assert!(decompressed.contains("@read1"));
        assert!(decompressed.contains("AAAA"));
        assert!(decompressed.contains("@read2"));
        assert!(decompressed.contains("GGGG"));
    }

    #[test]
    fn test_write_output_fastq_bzip() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("output.fastq");
        let (tx, rx) = crossbeam::channel::unbounded();
        let record1 = fastq::Record::new(
            fastq::record::Definition::new("read1", "read1"),
            "AAAA",
            "!!!!",
        );
        let record2 = fastq::Record::new(
            fastq::record::Definition::new("read2", "read2"),
            "GGGG",
            "!!!!",
        );
        tx.send(record1).unwrap();
        tx.send(record2).unwrap();
        drop(tx);
        let read_count = write_output_fastq(
            rx,
            &file_path,
            Some(niffler::compression::Format::Bzip),
            niffler::Level::One,
        )
        .unwrap();
        let reader =
            niffler::get_reader(Box::new(BufReader::new(File::open(&file_path).unwrap()))).unwrap();
        let mut decompressed = String::new();
        let mut decompressed_reader = BufReader::new(reader.0);
        decompressed_reader
            .read_to_string(&mut decompressed)
            .unwrap();

        assert_eq!(read_count, 2);
        assert!(decompressed.contains("@read1"));
        assert!(decompressed.contains("AAAA"));
        assert!(decompressed.contains("@read2"));
        assert!(decompressed.contains("GGGG"));
    }

    #[test]
    fn test_write_output_fasta_non_compressed() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("output.fasta");
        let (tx, rx) = crossbeam::channel::unbounded();
        let record1 = fastq::Record::new(
            fastq::record::Definition::new("read1", "read1"),
            "AAAA",
            "!!!!",
        );
        let record2 = fastq::Record::new(
            fastq::record::Definition::new("read2", "read2"),
            "GGGG",
            "!!!!",
        );
        tx.send(record1).unwrap();
        tx.send(record2).unwrap();
        drop(tx);
        let read_count = write_output_fasta(rx, &file_path).unwrap();
        let file_content = fs::read_to_string(file_path).unwrap();

        assert_eq!(read_count, 2);
        assert!(file_content.contains(">read1"));
        assert!(file_content.contains("AAAA"));
        assert!(file_content.contains(">read2"));
        assert!(file_content.contains("GGGG"));
    }

    #[test]
    fn test_write_output_fastq_error() {
        let file_path = PathBuf::from("/noperms.fastq");
        let (_, rx) = crossbeam::channel::unbounded();
        let result = write_output_fastq(
            rx,
            &file_path,
            Some(niffler::compression::Format::No),
            niffler::Level::One,
        );

        assert!(result.is_err());
    }

    #[test]
    fn test_write_output_fasta_file_creation_error() {
        let file_path = PathBuf::from("/noperms.fasta");
        let (_, rx) = crossbeam::channel::unbounded();
        let result = write_output_fasta(rx, &file_path);

        assert!(result.is_err());
    }

    #[test]
    fn test_output_subdirs_dont_exist() {
        let dir = tempdir().unwrap();
        let subdir = dir.path().join("subdir1");
        let subdir = subdir.join("subdir2");
        let file_path = subdir.join("output.fastq");
        assert!(!subdir.exists());
        let (tx, rx) = crossbeam::channel::unbounded();
        let record1 = fastq::Record::new(
            fastq::record::Definition::new("read1", "read1"),
            "AAAA",
            "!!!!",
        );
        tx.send(record1).unwrap();
        drop(tx);

        let read_count = write_output_fastq(
            rx,
            &file_path,
            Some(niffler::compression::Format::No),
            niffler::Level::One,
        )
        .unwrap();

        assert_eq!(read_count, 1);
        assert!(subdir.exists());
    }

    #[test]
    fn test_output_subdirs_already_exist() {
        let dir = tempdir().unwrap();
        let subdir = dir.path().join("subdir1");
        let subdir = subdir.join("subdir2");
        let file_path = subdir.join("output.fastq");
        fs::create_dir_all(&subdir).unwrap();
        assert!(subdir.exists());
        let (tx, rx) = crossbeam::channel::unbounded();
        let record1 = fastq::Record::new(
            fastq::record::Definition::new("read1", "read1"),
            "AAAA",
            "!!!!",
        );
        tx.send(record1).unwrap();
        drop(tx);

        let read_count = write_output_fastq(
            rx,
            &file_path,
            Some(niffler::compression::Format::No),
            niffler::Level::One,
        )
        .unwrap();

        assert_eq!(read_count, 1);
        assert!(subdir.exists());
    }
}
