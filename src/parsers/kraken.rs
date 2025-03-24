use anyhow::{anyhow, bail, Context, Result};
use fxhash::FxHashSet;
use log::{debug, info};
use std::collections::{HashMap, HashSet};
use std::fs;
use std::io::{BufRead, BufReader};

#[derive(Debug, Clone)]
pub struct Tree {
    pub taxon_id: i32,
    pub level_num: usize,
    pub children: Vec<usize>,
    pub parent: Option<usize>,
}

impl Tree {
    pub fn new(taxon_id: i32, level_num: usize, parent: Option<usize>) -> Tree {
        Tree {
            taxon_id,
            level_num,
            children: Vec::new(),
            parent,
        }
    }
}

#[derive(Debug, Clone)]
pub struct KrakenRecord {
    pub is_classified: bool,
    pub read_id: Vec<u8>,
    pub taxon_id: i32,
    pub length: String,
    pub lca_map: String,
}

#[derive(Debug, Clone)]
pub struct KrakenReportRecord {
    pub percent: f32,
    pub fragments_clade_rooted: i32,
    pub fragments_taxon: i32,
    pub rank: String,
    pub taxon_id: i32,
    pub level: usize,
    pub name: String,
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
/// A KrakenRecord containing the extracted information
fn process_kraken_output_line(kraken_output: &str) -> Result<KrakenRecord> {
    let mut fields = kraken_output.split('\t');

    let classification = fields
        .next()
        .ok_or_else(|| anyhow!("Missing classification field"))?;
    let read_id = fields
        .next()
        .ok_or_else(|| anyhow!("Missing read ID field"))?;
    let taxon_id = fields
        .next()
        .ok_or_else(|| anyhow!("Missing taxon ID field"))?;
    let length = fields
        .next()
        .ok_or_else(|| anyhow!("Missing length field"))?;
    let lca_map = fields.next().ok_or_else(|| anyhow!("Missing LCA map"))?;

    if fields.next().is_none() {
        let is_classified = classification == "C";
        let taxon_id = taxon_id
            .trim()
            .parse::<i32>()
            .with_context(|| format!("Error parsing taxon ID: '{}'", taxon_id))?;
        Ok(KrakenRecord {
            is_classified,
            read_id: read_id.as_bytes().to_vec(),
            taxon_id,
            length: length.to_string(),
            lca_map: lca_map.to_string(),
        })
    } else {
        bail!("Invalid kraken output line format: Expected 5 tab-separated fields, but got more");
    }
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
pub fn process_kraken_output(
    kraken_path: &str,
    exclude: bool,
    taxon_ids_to_save: &[i32],
) -> Result<FxHashSet<Vec<u8>>> {
    info!("Processing kraken output...");
    let taxon_ids_to_save: HashSet<i32> = taxon_ids_to_save.iter().cloned().collect();
    let mut reads_to_save = FxHashSet::default();
    let kraken_file = fs::File::open(kraken_path)
        .with_context(|| format!("Failed to open kraken output file: {}", kraken_path))?;
    let reader = BufReader::new(kraken_file);

    for line_result in reader.lines() {
        let line = line_result.context("Error reading kraken output line")?;
        let record = process_kraken_output_line(&line)?;
        if (exclude && !taxon_ids_to_save.contains(&record.taxon_id))
            || (!exclude && taxon_ids_to_save.contains(&record.taxon_id))
        {
            reads_to_save.insert(record.read_id);
        }
    }
    // let mut taxon_id_count = TAXON_ID_COUNT
    //     .lock()
    //     .map_err(|e| anyhow!("Failed to lock TAXON_ID_COUNT mutex: {}", e))?;
    // let mut taxon_ids = TAXON_IDS
    //     .lock()
    //     .map_err(|e| anyhow!("Failed to lock TAXON_IDS mutex: {}", e))?;
    // let mut total_read_count = TOTAL_READS
    //     .lock()
    //     .map_err(|e| anyhow!("Failed to lock TOTAL_READS mutex: {}", e))?;
    // let mut reads_to_extract = READS_TO_EXTRACT
    //     .lock()
    //     .map_err(|e| anyhow!("Failed to lock READS_TO_EXTRACT mutex: {}", e))?;
    //
    // *taxon_id_count = taxon_ids_to_save.len();
    // *taxon_ids = taxon_ids_to_save;
    // *total_read_count = total_reads; // Update this with the actual total_reads value
    // *reads_to_extract = reads_to_save.len();

    Ok(reads_to_save)
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
fn process_kraken_report_line(kraken_report: &str) -> Result<KrakenReportRecord> {
    let mut fields = kraken_report.split('\t');

    let percent_field = fields
        .next()
        .ok_or_else(|| anyhow!("Missing percent field"))?;
    let fragments_clade_rooted_field = fields
        .next()
        .ok_or_else(|| anyhow!("Missing fragments clade rooted field"))?;
    let fragments_taxon_field = fields
        .next()
        .ok_or_else(|| anyhow!("Missing fragments taxon field"))?;
    let rank_field = fields.next().ok_or_else(|| anyhow!("Missing rank field"))?;
    let taxon_id_field = fields
        .next()
        .ok_or_else(|| anyhow!("Missing taxon ID field"))?;
    let name_field = fields
        .next()
        .ok_or_else(|| anyhow!("Missing taxon name field"))?;
    if fields.next().is_none() {
        let percent = percent_field
            .trim()
            .parse::<f32>()
            .with_context(|| format!("Error parsing percent value: '{}'", percent_field))?;

        let fragments_clade_rooted = fragments_clade_rooted_field
            .trim()
            .parse::<i32>()
            .with_context(|| {
                format!(
                    "Error parsing fragments clade rooted: '{}'",
                    fragments_clade_rooted_field
                )
            })?;

        let fragments_taxon = fragments_taxon_field
            .trim()
            .parse::<i32>()
            .with_context(|| {
                format!("Error parsing fragments taxon: '{}'", fragments_taxon_field)
            })?;

        let taxon_id = taxon_id_field
            .trim()
            .parse::<i32>()
            .with_context(|| format!("Error parsing taxon ID: '{}'", taxon_id_field))?;

        let level = name_field.chars().take_while(|&c| c == ' ').count() / 2;

        Ok(KrakenReportRecord {
            percent,
            fragments_clade_rooted,
            fragments_taxon,
            rank: rank_field.to_string(),
            taxon_id,
            level,
            name: name_field.to_string(),
        })
    } else {
        bail!("Invalid kraken report line format: Expected 6 tab-separated fields, but got more");
    }
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
pub fn build_tree_from_kraken_report(
    taxon_to_save: i32,
    report_path: &str,
) -> Result<(Vec<Tree>, HashMap<i32, usize>)> {
    debug!("Building taxonomic tree from kraken report");
    // will store the tree
    let mut nodes = Vec::new();
    // taxonid -> index in the nodes vector
    let mut taxon_map = HashMap::new();

    let report_file = fs::File::open(report_path)
        .with_context(|| format!("Failed to open kraken report file: {}", report_path))?;

    let reader = BufReader::new(report_file);
    let mut prev_index = None;

    for line in reader.lines() {
        let line = line.context("Error reading kraken report line")?;
        let record = process_kraken_report_line(&line)?;
        // if taxon_id == 0, it's an unclassified read so we can skip
        if record.taxon_id == 0 {
            continue;
        }
        // 1 will be the root of the tree
        if record.taxon_id == 1 {
            let root_node = Tree::new(record.taxon_id, record.level, None);
            prev_index = Some(nodes.len());
            nodes.push(root_node);
        }
        // if the current level is not the same as the previous level + 1, then we are not at the correct parent, and need to move up the tree
        while let Some(parent_index) = prev_index {
            if record.level != nodes[parent_index].level_num + 1 {
                prev_index = nodes[parent_index].parent;
            } else {
                break;
            }
        }
        // once we have the correct parent, we can add the current node to the tree
        let curr_node = Tree::new(record.taxon_id, record.level, prev_index);
        let curr_index = nodes.len();
        nodes.push(curr_node);

        // add the current node
        if let Some(parent_index) = prev_index {
            nodes[parent_index].children.push(curr_index);
        }

        prev_index = Some(curr_index);

        // if the current taxon is one we want to save, add it to the map
        if record.taxon_id == taxon_to_save {
            taxon_map.insert(record.taxon_id, curr_index);
        }
    }

    if !taxon_map.contains_key(&taxon_to_save) {
        bail!("Taxon ID {} not found in the kraken report", taxon_to_save);
    }

    Ok((nodes, taxon_map))
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
pub fn extract_parents(
    taxon_map: &HashMap<i32, usize>,
    nodes: &[Tree],
    taxon_id: i32,
) -> Result<Vec<i32>> {
    // Backtracking traversal from the given taxon_id to the root
    let start_index = taxon_map
        .get(&taxon_id)
        .ok_or_else(|| anyhow!("Taxon ID {} not found in taxonomy map", taxon_id))?;

    let mut parents = Vec::new();
    parents.push(taxon_id);
    let mut curr_index = *start_index;

    while let Some(parent_index) = nodes[curr_index].parent {
        if parent_index >= nodes.len() {
            bail!(
                "Invalid parent index {} for node at index {}",
                parent_index,
                curr_index
            );
        }
        parents.push(nodes[parent_index].taxon_id);
        curr_index = parent_index;
    }

    Ok(parents)
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
pub fn extract_children(nodes: &[Tree], start_index: usize, result: &mut Vec<i32>) -> Result<()> {
    // recursive post-order traversal of the tree
    if start_index >= nodes.len() {
        bail!(
            "Invalid start index {} for node in tree of length {}",
            start_index,
            nodes.len()
        );
    }
    for &child_index in &nodes[start_index].children {
        if child_index >= nodes.len() {
            bail!(
                "Invalid child index {} for node at index {}",
                child_index,
                start_index
            );
        }
        extract_children(nodes, child_index, result)?;
    }
    result.push(nodes[start_index].taxon_id);
    Ok(())
}
