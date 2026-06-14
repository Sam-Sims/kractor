use color_eyre::{
    Result,
    eyre::{bail, ensure},
};
use fxhash::{FxHashMap, FxHashSet};
use log::info;
use serde::{Deserialize, Serialize};

use crate::{
    Cli,
    extract::{self, KractorResult, process_paired_end, process_single_end},
    parsers::{self, kraken::ProcessedKrakenOutput},
};

pub fn run(args: Cli) -> Result<()> {
    Kractor::new(args).run_inner()
}

struct Kractor {
    args: Cli,
    taxon_ids: Vec<i32>,
    missing_taxon_ids: Vec<i32>,
    reads_to_save: FxHashSet<Vec<u8>>,
    reads_per_taxon: FxHashMap<i32, usize>,
    summary: Option<Summary>,
}

impl Kractor {
    fn new(args: Cli) -> Self {
        Self {
            args,
            taxon_ids: Vec::new(),
            missing_taxon_ids: Vec::new(),
            reads_to_save: FxHashSet::default(),
            reads_per_taxon: FxHashMap::default(),
            summary: None,
        }
    }

    fn run_inner(&mut self) -> Result<()> {
        info!(
            "Starting kractor at {}",
            chrono::Local::now().format("%H:%M:%S")
        );
        self.validate_outputs()?;
        self.collect_taxa()?;
        info!("Processing Kraken2 output file");
        self.process_kraken_output()?;
        info!("Processing reads");
        self.process_reads()?;
        info!("Complete at {}", chrono::Local::now().format("%H:%M:%S"));
        self.output_summary()?;
        Ok(())
    }

    fn validate_outputs(&self) -> Result<()> {
        for out_file in &self.args.output {
            ensure!(
                !out_file.exists(),
                "Output file already exists: {}",
                out_file.display()
            );
        }
        Ok(())
    }

    fn collect_taxa(&mut self) -> Result<()> {
        let collected = extract::collect_taxa_to_save(
            self.args.report.as_deref(),
            self.args.children,
            self.args.parents,
            &self.args.taxid,
            !self.args.no_report_header_detect,
        )?;
        self.taxon_ids = collected.found;
        self.missing_taxon_ids = collected.missing;
        Ok(())
    }

    fn process_kraken_output(&mut self) -> Result<()> {
        let ProcessedKrakenOutput {
            reads_to_save,
            reads_per_taxon,
        } = parsers::kraken::process_kraken_output(
            &self.args.kraken,
            self.args.exclude,
            &self.taxon_ids,
        )?;
        self.reads_to_save = reads_to_save;
        self.reads_per_taxon = reads_per_taxon;

        if self.reads_to_save.is_empty() {
            bail!("No reads found for the specified taxon ID(s). Nothing to extract.");
        }

        info!("Identified {} reads to save", self.reads_to_save.len());
        Ok(())
    }

    fn process_reads(&mut self) -> Result<()> {
        let paired = self.args.input.len() == 2;
        let input_layout = if paired { "paired" } else { "single" };
        let reads_extracted_per_taxon = self.get_reads_extracted_per_taxon();

        let result = if paired {
            let (r1, r2) = process_paired_end(
                &self.reads_to_save,
                &self.args.input,
                &self.args.output,
                self.args.output_type,
                self.args.compression_level,
                self.args.output_format,
            )?;

            KractorResult {
                reads_parsed: r1.reads_parsed + r2.reads_parsed,
                reads_output: r1.reads_output + r2.reads_output,
                input_format: r1.input_format,
                output_format: r1.output_format,
            }
        } else {
            process_single_end(
                &self.reads_to_save,
                &self.args.input,
                &self.args.output,
                self.args.output_type,
                self.args.compression_level,
                self.args.output_format,
            )?
        };

        self.summary = Some(Summary {
            kractor_version: env!("CARGO_PKG_VERSION").to_string(),
            input_layout: input_layout.to_string(),
            input_sequence_format: result.input_format.to_string(),
            output_sequence_format: result.output_format.to_string(),
            requested_taxon_ids: self.args.taxid.clone(),
            matched_taxon_ids: self.taxon_ids.clone(),
            requested_taxon_ids_not_found: self.missing_taxon_ids.clone(),
            total_input_records: result.reads_parsed,
            total_output_records: result.reads_output,
            extraction_fraction: result.reads_output as f64 / result.reads_parsed as f64,
            assigned_reads_per_taxon: reads_extracted_per_taxon,
        });

        Ok(())
    }

    fn output_summary(&self) -> Result<()> {
        if self.args.summary
            && let Some(summary) = &self.summary
        {
            let json = serde_json::to_string_pretty(summary)?;
            println!("{json}");
        }
        Ok(())
    }

    fn get_reads_extracted_per_taxon(&self) -> FxHashMap<i32, usize> {
        let mut reads_extracted_per_taxon = self.reads_per_taxon.clone();
        for taxon_id in &self.taxon_ids {
            reads_extracted_per_taxon.entry(*taxon_id).or_insert(0);
        }
        reads_extracted_per_taxon
    }
}

#[derive(Serialize, Deserialize)]
struct Summary {
    kractor_version: String,
    input_layout: String,
    input_sequence_format: String,
    output_sequence_format: String,
    requested_taxon_ids: Vec<i32>,
    matched_taxon_ids: Vec<i32>,
    requested_taxon_ids_not_found: Vec<i32>,
    total_input_records: usize,
    total_output_records: usize,
    extraction_fraction: f64,
    assigned_reads_per_taxon: FxHashMap<i32, usize>,
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use tempfile::tempdir;

    use super::*;
    use crate::cli::OutputFormat;

    #[test]
    fn test_output_doesnt_exist() {
        let temp_dir = tempdir().unwrap();
        let output_file = temp_dir.path().join("output.fastq");
        let input_files = vec![PathBuf::from("input1.fastq"), PathBuf::from("input2.fastq")];
        let args = Cli {
            input: input_files,
            output: vec![output_file],
            kraken: PathBuf::from("kraken_output.txt"),
            report: None,
            taxid: vec![1, 2, 3],
            output_type: None,
            compression_level: niffler::Level::One,
            parents: false,
            children: false,
            exclude: false,
            output_format: OutputFormat::Auto,
            summary: false,
            no_report_header_detect: false,
            verbose: false,
        };
        let kractor = Kractor::new(args);
        assert!(kractor.validate_outputs().is_ok());
    }

    #[test]
    fn test_output_exists() {
        let temp_dir = tempdir().unwrap();
        let output_file = temp_dir.path().join("output.fastq");
        std::fs::File::create(&output_file).unwrap();
        let input_files = vec![PathBuf::from("input.fastq")];
        let args = Cli {
            input: input_files,
            output: vec![output_file],
            kraken: PathBuf::from("kraken_output.txt"),
            report: None,
            taxid: vec![1, 2, 3],
            output_type: None,
            compression_level: niffler::Level::One,
            parents: false,
            children: false,
            exclude: false,
            output_format: OutputFormat::Auto,
            summary: false,
            no_report_header_detect: false,
            verbose: false,
        };
        let kractor = Kractor::new(args);
        assert!(kractor.validate_outputs().is_err());
    }

    #[test]
    fn test_get_reads_extracted_per_taxon() {
        let input_files = vec![PathBuf::from("input.fastq")];
        let args = Cli {
            input: input_files,
            output: vec![PathBuf::from("output.fastq")],
            kraken: PathBuf::from("kraken_output.txt"),
            report: None,
            taxid: vec![2901879, 227984],
            output_type: None,
            compression_level: niffler::Level::One,
            parents: false,
            children: false,
            exclude: false,
            output_format: OutputFormat::Auto,
            summary: false,
            no_report_header_detect: false,
            verbose: false,
        };
        let mut kractor = Kractor::new(args);
        kractor.taxon_ids = vec![2901879, 227984];
        kractor.reads_per_taxon.insert(227984, 257);

        let reads_extracted_per_taxon = kractor.get_reads_extracted_per_taxon();

        assert_eq!(reads_extracted_per_taxon.get(&2901879), Some(&0));
        assert_eq!(reads_extracted_per_taxon.get(&227984), Some(&257));
    }
}
