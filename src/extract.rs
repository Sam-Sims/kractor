use crate::parsers::fastx::{parse_fastq, write_output_fasta, write_output_fastq};
use crate::parsers::kraken::{build_tree_from_kraken_report, extract_children, extract_parents};
use color_eyre::{eyre::bail, eyre::eyre, eyre::WrapErr, Result};
use crossbeam::{channel, thread};
use fxhash::FxHashSet;
use log::debug;
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
    input: &[PathBuf],
    output: &[PathBuf],
    compression_type: Option<niffler::Format>,
    compression_level: niffler::Level,
    fasta: bool,
) -> Result<usize> {
    thread::scope(|scope| -> Result<usize> {
        let (tx, rx) = channel::unbounded::<fastq::Record>();

        let reader = scope.spawn(|_| {
            let result = parse_fastq(&input[0], reads_to_save, &tx);
            drop(tx);
            result.wrap_err_with(|| format!("Failed to parse input file: {:?}", input[0]))
        });

        let writer = scope.spawn(|_| {
            if !fasta {
                write_output_fastq(rx, &output[0], compression_type, compression_level)
                    .wrap_err_with(|| format!("Failed to write output file: {:?}", output[0]))
            } else {
                write_output_fasta(rx, &output[0])
                    .wrap_err_with(|| format!("Failed to write output file: {:?}", output[0]))
            }
        });

        reader
            .join()
            .map_err(|_| eyre!("Reader thread panicked"))??;
        let total_reads_output = writer
            .join()
            .map_err(|_| eyre!("Writer thread panicked"))??;
        Ok(total_reads_output)
    })
    .map_err(|_| eyre!("Thread communication error"))?
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
    input: &[PathBuf],
    output: &[PathBuf],
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
            result.wrap_err_with(|| format!("Failed to parse first input file: {:?}", input[0]))
        });

        let reader2 = scope.spawn(|_| {
            let result = parse_fastq(&input[1], reads_to_save, &tx2);
            drop(tx2);
            result.wrap_err_with(|| format!("Failed to parse second input file: {:?}", input[1]))
        });

        let writer1 = scope.spawn(|_| {
            if !fasta {
                write_output_fastq(rx1, &output[0], compression_type, compression_level)
                    .wrap_err_with(|| {
                        format!(
                            "Failed to write FASTQ output to first file: {:?}",
                            output[0]
                        )
                    })
            } else {
                write_output_fasta(rx1, &output[0]).wrap_err_with(|| {
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
                    .wrap_err_with(|| {
                        format!(
                            "Failed to write FASTQ output to second file: {:?}",
                            output[1]
                        )
                    })
            } else {
                write_output_fasta(rx2, &output[1]).wrap_err_with(|| {
                    format!(
                        "Failed to write FASTA output to second file: {:?}",
                        output[1]
                    )
                })
            }
        });

        reader1
            .join()
            .map_err(|_| eyre!("Reader thread for file1 panicked"))??;
        reader2
            .join()
            .map_err(|_| eyre!("Reader thread for file2 panicked"))??;
        let total_reads_output_pair1 = writer1
            .join()
            .map_err(|_| eyre!("Writer thread for file1 panicked"))??;
        let total_reads_output_pair2 = writer2
            .join()
            .map_err(|_| eyre!("Writer thread for file2 panicked"))??;
        Ok((total_reads_output_pair1, total_reads_output_pair2))
    })
    .map_err(|_| eyre!("Thread communication error"))?
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
        return Err(eyre!("Report required when parents or children is enabled"));
    }

    if let Some(report_path) = report {
        let (nodes, taxon_map) = build_tree_from_kraken_report(&taxids, report_path)
            .wrap_err("Failed to build tree from Kraken report")?;

        if children {
            debug!("Extracting children");
            let mut children = Vec::new();
            for taxid in &taxids {
                if let Some(&node_index) = taxon_map.get(taxid) {
                    extract_children(&nodes, node_index, &mut children).wrap_err_with(|| {
                        format!("Failed to extract children for taxon ID {}", taxid)
                    })?;
                } else {
                    return Err(eyre!("Taxon ID {} not found in taxonomy map", taxid));
                }
            }
            taxon_ids_to_save.extend(&children);
        } else if parents {
            debug!("Extracting parents");
            for taxid in &taxids {
                taxon_ids_to_save.extend(
                    extract_parents(&taxon_map, &nodes, *taxid).wrap_err_with(|| {
                        format!("Failed to extract parents for taxon ID {}", taxid)
                    })?,
                );
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

#[cfg(test)]

mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write;
    use tempfile::tempdir;

    #[test]
    fn test_process_single_end_fastq() {
        let dir = tempdir().unwrap();
        let input_path = dir.path().join("input.fastq");
        let output_path = dir.path().join("output.fastq");
        let test_data = "@read1\nAAAA\n+\n!!!!\n@read2\nGGGG\n+\n!!!!\n";
        let mut file = File::create(&input_path).unwrap();
        file.write_all(test_data.as_bytes()).unwrap();
        let mut reads_to_save = FxHashSet::default();
        reads_to_save.insert(b"read1".to_vec());
        let input = vec![input_path];
        let output = vec![output_path.clone()];
        let read_count = process_single_end(
            &reads_to_save,
            &input,
            &output,
            Some(niffler::compression::Format::No),
            niffler::Level::One,
            false,
        )
        .unwrap();
        let file_content = std::fs::read_to_string(output_path).unwrap();

        assert_eq!(read_count, 1);
        assert!(file_content.contains("@read1"));
        assert!(file_content.contains("AAAA"));
        assert!(!file_content.contains("@read2"));
    }

    #[test]
    fn test_process_single_end_fasta() {
        let dir = tempdir().unwrap();
        let input_path = dir.path().join("input.fastq");
        let output_path = dir.path().join("output.fastq");
        let test_data = "@read1\nAAAA\n+\n!!!!\n@read2\nGGGG\n+\n!!!!\n";
        let mut file = File::create(&input_path).unwrap();
        file.write_all(test_data.as_bytes()).unwrap();
        let mut reads_to_save = FxHashSet::default();
        reads_to_save.insert(b"read1".to_vec());
        let input = vec![input_path];
        let output = vec![output_path.clone()];
        let read_count = process_single_end(
            &reads_to_save,
            &input,
            &output,
            Some(niffler::compression::Format::No),
            niffler::Level::One,
            true,
        )
        .unwrap();
        let file_content = std::fs::read_to_string(output_path).unwrap();

        assert_eq!(read_count, 1);
        assert!(file_content.contains(">read1"));
        assert!(file_content.contains("AAAA"));
        assert!(!file_content.contains("@read2"));
    }

    #[test]
    fn test_process_single_end_not_found() {
        let nonexistent_path = PathBuf::from("idontexist.fastq");
        let output_path = PathBuf::from("output.fastq");
        let reads_to_save = FxHashSet::default();
        let input = vec![nonexistent_path];
        let output = vec![output_path];

        let result = process_single_end(
            &reads_to_save,
            &input,
            &output,
            None,
            niffler::Level::One,
            false,
        );

        assert!(result.is_err());
    }

    #[test]
    fn test_process_paired_end_fastq() {
        let dir = tempdir().unwrap();
        let input_path1 = dir.path().join("input1.fastq");
        let input_path2 = dir.path().join("input2.fastq");
        let output_path1 = dir.path().join("output1.fastq");
        let output_path2 = dir.path().join("output2.fastq");
        let test_data1 = "@read1\nAAAA\n+\n!!!!\n@read2\nGGGG\n+\n!!!!\n";
        let test_data2 = "@read1\nTTTT\n+\n!!!!\n@read2\nCCCC\n+\n!!!!\n";
        let mut file1 = File::create(&input_path1).unwrap();
        file1.write_all(test_data1.as_bytes()).unwrap();
        let mut file2 = File::create(&input_path2).unwrap();
        file2.write_all(test_data2.as_bytes()).unwrap();
        let mut reads_to_save = FxHashSet::default();
        reads_to_save.insert(b"read1".to_vec());
        let input = vec![input_path1, input_path2];
        let output = vec![output_path1.clone(), output_path2.clone()];
        let (read_count1, read_count2) = process_paired_end(
            &reads_to_save,
            &input,
            &output,
            Some(niffler::compression::Format::No),
            niffler::Level::One,
            false,
        )
        .unwrap();
        let file_content1 = std::fs::read_to_string(output_path1).unwrap();
        let file_content2 = std::fs::read_to_string(output_path2).unwrap();

        assert_eq!(read_count1, 1);
        assert_eq!(read_count2, 1);
        assert!(file_content1.contains("@read1"));
        assert!(file_content1.contains("AAAA"));
        assert!(!file_content1.contains("@read2"));
        assert!(file_content2.contains("@read1"));
        assert!(file_content2.contains("TTTT"));
    }

    #[test]
    fn test_process_paired_end_fasta() {
        let dir = tempdir().unwrap();
        let input_path1 = dir.path().join("input1.fastq");
        let input_path2 = dir.path().join("input2.fastq");
        let output_path1 = dir.path().join("output1.fasta");
        let output_path2 = dir.path().join("output2.fasta");
        let test_data1 = "@read1\nAAAA\n+\n!!!!\n@read2\nGGGG\n+\n!!!!\n";
        let test_data2 = "@read1\nTTTT\n+\n!!!!\n@read2\nCCCC\n+\n!!!!\n";
        let mut file1 = File::create(&input_path1).unwrap();
        file1.write_all(test_data1.as_bytes()).unwrap();
        let mut file2 = File::create(&input_path2).unwrap();
        file2.write_all(test_data2.as_bytes()).unwrap();
        let mut reads_to_save = FxHashSet::default();
        reads_to_save.insert(b"read1".to_vec());
        let input = vec![input_path1, input_path2];
        let output = vec![output_path1.clone(), output_path2.clone()];
        let (read_count1, read_count2) = process_paired_end(
            &reads_to_save,
            &input,
            &output,
            Some(niffler::compression::Format::No),
            niffler::Level::One,
            true,
        )
        .unwrap();
        let file_content1 = std::fs::read_to_string(output_path1).unwrap();
        let file_content2 = std::fs::read_to_string(output_path2).unwrap();

        assert_eq!(read_count1, 1);
        assert_eq!(read_count2, 1);
        assert!(file_content1.contains(">read1"));
        assert!(file_content1.contains("AAAA"));
        assert!(!file_content1.contains("@read2"));
        assert!(file_content2.contains(">read1"));
        assert!(file_content2.contains("TTTT"));
    }

    fn create_test_kraken_report(dir: &tempfile::TempDir) -> PathBuf {
        let report_path = dir.path().join("report.txt");
        let test_data = "\
        21.36\t745591\t745591\tU\t0\tunclassified
        78.64\t2745487\t1646\tR\t1\troot
        78.58\t2743340\t1360\tR1\t131567\t  cellular organisms
        78.21\t2730479\t8458\tD\t2\t    Bacteria
        61.55\t2148918\t1359\tD1\t1783272\t      Terrabacteria group
        61.40\t2143487\t321\tP\t1239\t        Bacillota
        61.37\t2142480\t8314\tC\t91062\t          Bacilli2
        61.37\t2142480\t8314\tC\t91061\t          Bacilli
        38.95\t1359681\t1300\tO\t1385\t            Bacillales
        16.53\t577203\t366\tF\t186817\t              Bacillaceae
        16.50\t576156\t22486\tG\t1386\t                Bacillus";

        let mut file = File::create(&report_path).unwrap();
        file.write_all(test_data.as_bytes()).unwrap();
        report_path
    }

    #[test]
    fn test_error_when_no_report_and_parents_or_children() {
        let result = collect_taxons_to_save(&None, true, false, vec![1]);
        assert!(result.is_err());
        let result = collect_taxons_to_save(&None, false, true, vec![1]);
        assert!(result.is_err());
    }

    #[test]
    fn test_no_report() {
        let taxids = vec![123, 456, 789];
        let saved_taxons = collect_taxons_to_save(&None, false, false, taxids.clone()).unwrap();

        assert_eq!(saved_taxons, taxids);
    }

    #[test]
    fn test_no_parents_no_children() {
        let dir = tempdir().unwrap();
        let report_path = create_test_kraken_report(&dir);
        let taxids = vec![1385, 1386, 91061];
        let saved_taxons =
            collect_taxons_to_save(&Some(report_path), false, false, taxids.clone()).unwrap();

        assert_eq!(saved_taxons, taxids);
    }

    #[test]
    fn test_children() {
        let dir = tempdir().unwrap();
        let report_path = create_test_kraken_report(&dir);
        let taxids = vec![1239];
        let saved_taxons = collect_taxons_to_save(&Some(report_path), true, false, taxids).unwrap();

        assert!(saved_taxons.contains(&1239));
        assert!(saved_taxons.contains(&91062));
        assert!(saved_taxons.contains(&91061));
    }

    #[test]
    fn test_parents() {
        let dir = tempdir().unwrap();
        let report_path = create_test_kraken_report(&dir);
        let taxids = vec![91061];
        let saved_taxons = collect_taxons_to_save(&Some(report_path), false, true, taxids).unwrap();

        assert!(saved_taxons.contains(&91061));
        assert!(saved_taxons.contains(&1239));
        assert!(saved_taxons.contains(&1783272));
        assert!(saved_taxons.contains(&131567));
        assert!(saved_taxons.contains(&2));
    }

    #[test]
    fn test_taxon_not_exist() {
        let dir = tempdir().unwrap();
        let report_path = create_test_kraken_report(&dir);
        let taxids = vec![999];
        let result = collect_taxons_to_save(&Some(report_path), true, false, taxids);

        assert!(result.is_err());
    }

    #[test]
    fn test_dedup_and_sort() {
        let taxids = vec![456, 123, 456, 789, 123];
        let saved_taxons = collect_taxons_to_save(&None, false, false, taxids).unwrap();

        assert_eq!(saved_taxons, vec![123, 456, 789]);
    }

    #[test]
    fn test_empty_result() {
        let result = collect_taxons_to_save(&None, false, false, vec![]);

        assert!(result.is_err());
    }
}
