use crate::extract::{process_paired_end, process_single_end};
use crate::{extract, parsers, Cli};
use color_eyre::Result;
use fxhash::FxHashSet;
use log::{debug, info};
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
struct Summary {
    taxon_count: usize,
    taxon_ids: Vec<i32>,
    reads_output: ReadCounts,
    input_format: String,
}

#[derive(Serialize, Deserialize)]
enum ReadCounts {
    Single {
        total: usize,
    },
    Paired {
        total: usize,
        read1: usize,
        read2: usize,
    },
}

pub struct Kractor {
    args: Cli,
    taxon_ids: Vec<i32>,
    reads_to_save: FxHashSet<Vec<u8>>,
    summary: Option<Summary>,
}

impl Kractor {
    pub fn new(args: Cli) -> Self {
        Self {
            args,
            taxon_ids: Vec::new(),
            reads_to_save: FxHashSet::default(),
            summary: None,
        }
    }

    fn collect_taxons(&mut self) -> Result<()> {
        self.taxon_ids = extract::collect_taxons_to_save(
            &self.args.report,
            self.args.children,
            self.args.parents,
            self.args.taxid.clone(),
        )?;
        debug!("Taxon IDs identified: {:?}", self.taxon_ids);
        Ok(())
    }

    fn process_kraken_output(&mut self) -> Result<()> {
        self.reads_to_save = parsers::kraken::process_kraken_output(
            &self.args.kraken,
            self.args.exclude,
            &self.taxon_ids,
        )?;

        debug!("Identified {} reads to save", self.reads_to_save.len());
        Ok(())
    }

    fn process_reads(&mut self) -> Result<()> {
        let paired = self.args.input.len() == 2;
        let input_format = if paired { "paired" } else { "single" };

        let reads_output = if paired {
            let (count1, count2) = process_paired_end(
                &self.reads_to_save,
                &self.args.input,
                &self.args.output,
                self.args.output_type,
                self.args.compression_level,
                self.args.output_fasta,
            )?;

            ReadCounts::Paired {
                total: count1 + count2,
                read1: count1,
                read2: count2,
            }
        } else {
            let total = process_single_end(
                &self.reads_to_save,
                &self.args.input,
                &self.args.output,
                self.args.output_type,
                self.args.compression_level,
                self.args.output_fasta,
            )?;

            ReadCounts::Single { total }
        };

        self.summary = Some(Summary {
            taxon_count: self.taxon_ids.len(),
            taxon_ids: self.taxon_ids.clone(),
            reads_output,
            input_format: input_format.to_string(),
        });

        Ok(())
    }

    fn output_summary(&self) -> Result<()> {
        if let Some(summary) = &self.summary {
            if self.args.json {
                let json = serde_json::to_string_pretty(summary)?;
                println!("{}", json);
            }
        }
        Ok(())
    }

    pub fn run(&mut self) -> Result<()> {
        info!(
            "Starting kractor at {}",
            chrono::Local::now().format("%H:%M:%S")
        );
        self.collect_taxons()?;
        info!("{} taxons identified to save", self.taxon_ids.len());
        info!("Processing Kraken2 output file");
        self.process_kraken_output()?;
        info!("Processing reads");
        self.process_reads()?;
        self.output_summary()?;
        info!("Complete at {}", chrono::Local::now().format("%H:%M:%S"));
        Ok(())
    }
}
