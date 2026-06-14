use crate::cli::OutputFormat;
use color_eyre::eyre::{Context, Result, eyre};
use crossbeam::channel::{Receiver, Sender};
use fxhash::FxHashSet;
use log::{debug, trace};
use std::path::Path;
use std::time::{Duration, Instant};
use std::{fs, io};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FastxFormat {
    Fasta,
    Fastq,
}

impl std::fmt::Display for FastxFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FastxFormat::Fasta => f.write_str("fasta"),
            FastxFormat::Fastq => f.write_str("fastq"),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastxRecord {
    pub id: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Option<Vec<u8>>,
}

pub fn parse_fastx(
    file_path: &Path,
    reads_to_save: &FxHashSet<Vec<u8>>,
    tx: &Sender<FastxRecord>,
) -> Result<(usize, FastxFormat)> {
    const PROGRESS_UPDATE_INTERVAL: Duration = Duration::from_millis(1500);

    let mut num_reads = 0;
    let mut input_format = None;
    let mut last_progress_update = Instant::now();

    let mut fastx_reader = needletail::parse_fastx_file(file_path)
        .wrap_err_with(|| format!("Failed to parse FASTX file: {}", file_path.display()))?;

    while let Some(result) = fastx_reader.next() {
        let record = result
            .wrap_err_with(|| format!("Error reading FASTX record at position {num_reads}"))?;

        let record_format = match record.format() {
            needletail::parser::Format::Fasta => FastxFormat::Fasta,
            needletail::parser::Format::Fastq => FastxFormat::Fastq,
        };
        input_format.get_or_insert(record_format);

        let record_id = record.id();
        let read_id = read_id(record_id);
        if reads_to_save.contains(&read_id.to_vec()) {
            tx.send(FastxRecord {
                id: record_id.to_vec(),
                seq: record.seq().into_owned(),
                qual: record.qual().map(Vec::from),
            })
            .wrap_err("Error sending record")?;
        }

        num_reads += 1;

        if last_progress_update.elapsed() >= PROGRESS_UPDATE_INTERVAL {
            trace!("Processed {num_reads} reads");
            last_progress_update = Instant::now();
        }
    }

    let input_format = input_format.ok_or_else(|| {
        eyre!(
            "No FASTA or FASTQ records found in input file: {}",
            file_path.display()
        )
    })?;

    Ok((num_reads, input_format))
}

pub fn detect_fastx_format(file_path: &Path) -> Result<FastxFormat> {
    let mut fastx_reader = needletail::parse_fastx_file(file_path)
        .wrap_err_with(|| format!("Failed to parse FASTX file: {}", file_path.display()))?;

    let record = fastx_reader
        .next()
        .ok_or_else(|| {
            eyre!(
                "No FASTA or FASTQ records found in input file: {}",
                file_path.display()
            )
        })?
        .wrap_err_with(|| {
            format!(
                "Error reading first FASTX record from {}",
                file_path.display()
            )
        })?;

    match record.format() {
        needletail::parser::Format::Fasta => Ok(FastxFormat::Fasta),
        needletail::parser::Format::Fastq => Ok(FastxFormat::Fastq),
    }
}

pub fn resolve_output_format(input: FastxFormat, requested: OutputFormat) -> FastxFormat {
    match requested {
        OutputFormat::Auto => input,
        OutputFormat::Fasta => FastxFormat::Fasta,
        OutputFormat::Fastq => FastxFormat::Fastq,
    }
}

pub fn write_output_fastx(
    rx: Receiver<FastxRecord>,
    out_file: &Path,
    output_format: FastxFormat,
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

    if let Some(parent) = out_file.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)
            .wrap_err_with(|| format!("Failed to create output directory: {}", parent.display()))?;
    }

    let out_file_handle = fs::File::create(out_file)
        .wrap_err_with(|| format!("Failed to create output file: {}", out_file.display()))?;

    let file_handle = Box::new(io::BufWriter::new(out_file_handle));
    let mut writer = niffler::get_writer(file_handle, compression_type, compression_level)
        .wrap_err("Failed to create niffler writer")?;

    for record in rx {
        match output_format {
            FastxFormat::Fasta => needletail::parser::write_fasta(
                &record.id,
                &record.seq,
                writer.as_mut(),
                needletail::parser::LineEnding::Unix,
            )
            .wrap_err_with(|| format!("Error writing FASTA record: {record:?}"))?,
            FastxFormat::Fastq => needletail::parser::write_fastq(
                &record.id,
                &record.seq,
                record.qual.as_deref(),
                writer.as_mut(),
                needletail::parser::LineEnding::Unix,
            )
            .wrap_err_with(|| format!("Error writing FASTQ record: {record:?}"))?,
        }

        read_output_count += 1;
    }

    Ok(read_output_count)
}

fn read_id(record_id: &[u8]) -> &[u8] {
    record_id
        .split(|byte| byte.is_ascii_whitespace())
        .next()
        .unwrap_or(record_id)
}

fn infer_compression(file_path: &Path) -> niffler::compression::Format {
    match file_path
        .extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| ext.to_ascii_lowercase())
    {
        Some(ref ext) if ext == "gz" => niffler::compression::Format::Gzip,
        Some(ref ext) if ext == "bz2" => niffler::compression::Format::Bzip,
        _ => niffler::compression::Format::No,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::{BufReader, Read, Write};
    use std::path::PathBuf;
    use tempfile::tempdir;

    fn fastx_record(id: &str, seq: &str, qual: Option<&str>) -> FastxRecord {
        FastxRecord {
            id: id.as_bytes().to_vec(),
            seq: seq.as_bytes().to_vec(),
            qual: qual.map(|value| value.as_bytes().to_vec()),
        }
    }

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
    fn test_infer_compression_no_extension() {
        let file_path = PathBuf::from("test");
        let compression = infer_compression(&file_path);

        assert_eq!(compression, niffler::compression::Format::No);
    }

    #[test]
    fn test_infer_compression_uppercase_extension() {
        let file_path = PathBuf::from("test.GZ");
        let compression = infer_compression(&file_path);

        assert_eq!(compression, niffler::compression::Format::Gzip);
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
        let (read_count, input_format) = parse_fastx(&file_path, &reads_to_save, &tx).unwrap();
        drop(tx);
        let results: Vec<FastxRecord> = rx.iter().collect();

        assert_eq!(read_count, 3);
        assert_eq!(input_format, FastxFormat::Fastq);
        assert_eq!(results.len(), 2);
        assert_eq!(results[0].id, b"read1");
        assert_eq!(results[1].id, b"read3");
        assert_eq!(results[0].seq, b"AAAA");
        assert_eq!(results[1].seq, b"TTTT");
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
        let (read_count, input_format) = parse_fastx(&file_path, &reads_to_save, &tx).unwrap();
        drop(tx);
        let results: Vec<FastxRecord> = rx.iter().collect();

        assert_eq!(read_count, 3);
        assert_eq!(input_format, FastxFormat::Fastq);
        assert_eq!(results.len(), 0);
    }

    #[test]
    fn test_parse_fastq_file_not_found() {
        let file_path = PathBuf::from("idontexist.fastq");
        let reads_to_save = FxHashSet::default();
        let (tx, _rx) = crossbeam::channel::unbounded();
        let result = parse_fastx(&file_path, &reads_to_save, &tx);

        assert!(result.is_err());
    }

    #[test]
    fn test_write_output_fastq_non_compressed() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("output.fastq");
        let (tx, rx) = crossbeam::channel::unbounded();
        tx.send(fastx_record("read1", "AAAA", Some("!!!!")))
            .unwrap();
        tx.send(fastx_record("read2", "GGGG", Some("!!!!")))
            .unwrap();
        drop(tx);
        let read_count = write_output_fastx(
            rx,
            &file_path,
            FastxFormat::Fastq,
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
        tx.send(fastx_record("read1", "AAAA", Some("!!!!")))
            .unwrap();
        tx.send(fastx_record("read2", "GGGG", Some("!!!!")))
            .unwrap();
        drop(tx);
        let read_count = write_output_fastx(
            rx,
            &file_path,
            FastxFormat::Fastq,
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
        tx.send(fastx_record("read1", "AAAA", Some("!!!!")))
            .unwrap();
        tx.send(fastx_record("read2", "GGGG", Some("!!!!")))
            .unwrap();
        drop(tx);
        let read_count = write_output_fastx(
            rx,
            &file_path,
            FastxFormat::Fastq,
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
        tx.send(fastx_record("read1", "AAAA", Some("!!!!")))
            .unwrap();
        tx.send(fastx_record("read2", "GGGG", Some("!!!!")))
            .unwrap();
        drop(tx);
        let read_count = write_output_fastx(
            rx,
            &file_path,
            FastxFormat::Fasta,
            Some(niffler::compression::Format::No),
            niffler::Level::One,
        )
        .unwrap();
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
        let result = write_output_fastx(
            rx,
            &file_path,
            FastxFormat::Fastq,
            Some(niffler::compression::Format::No),
            niffler::Level::One,
        );

        assert!(result.is_err());
    }

    #[test]
    fn test_write_output_fasta_file_creation_error() {
        let file_path = PathBuf::from("/noperms.fasta");
        let (_, rx) = crossbeam::channel::unbounded();
        let result = write_output_fastx(
            rx,
            &file_path,
            FastxFormat::Fasta,
            Some(niffler::compression::Format::No),
            niffler::Level::One,
        );

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
        tx.send(fastx_record("read1", "AAAA", Some("!!!!")))
            .unwrap();
        drop(tx);

        let read_count = write_output_fastx(
            rx,
            &file_path,
            FastxFormat::Fastq,
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
        tx.send(fastx_record("read1", "AAAA", Some("!!!!")))
            .unwrap();
        drop(tx);

        let read_count = write_output_fastx(
            rx,
            &file_path,
            FastxFormat::Fastq,
            Some(niffler::compression::Format::No),
            niffler::Level::One,
        )
        .unwrap();

        assert_eq!(read_count, 1);
        assert!(subdir.exists());
    }
}
