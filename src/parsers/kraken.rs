use color_eyre::{eyre::bail, eyre::eyre, eyre::Context, Result};
use fxhash::FxHashSet;
use log::debug;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

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

fn process_kraken_output_line(kraken_output: &str) -> Result<KrakenRecord> {
    let mut fields = kraken_output.split('\t');

    let classification = fields
        .next()
        .ok_or_else(|| eyre!("Missing classification field in the kraken output file"))?;
    let read_id = fields
        .next()
        .ok_or_else(|| eyre!("Missing read ID field in the kraken output file"))?;
    let taxon_id = fields
        .next()
        .ok_or_else(|| eyre!("Missing taxon ID field"))?;
    let length = fields
        .next()
        .ok_or_else(|| eyre!("Missing length field in the kraken output file"))?;
    let lca_map = fields
        .next()
        .ok_or_else(|| eyre!("Missing LCA map in the kraken output file"))?;

    if fields.next().is_none() {
        let is_classified = classification == "C";
        let taxon_id = taxon_id
            .trim()
            .parse::<i32>()
            .wrap_err_with(|| format!("Error parsing taxon ID: '{}'", taxon_id))?;
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

pub fn process_kraken_output(
    kraken_path: &PathBuf,
    exclude: bool,
    taxon_ids_to_save: &[i32],
) -> Result<FxHashSet<Vec<u8>>> {
    let taxon_ids_to_save: HashSet<i32> = taxon_ids_to_save.iter().cloned().collect();
    let mut reads_to_save = FxHashSet::default();
    let kraken_file = fs::File::open(kraken_path)
        .wrap_err_with(|| format!("Failed to open kraken output file: {:?}", kraken_path))?;
    let reader = BufReader::new(kraken_file);

    for line_result in reader.lines() {
        let line = line_result.wrap_err("Error reading kraken output line")?;
        let record = process_kraken_output_line(&line)?;
        if (exclude && !taxon_ids_to_save.contains(&record.taxon_id))
            || (!exclude && taxon_ids_to_save.contains(&record.taxon_id))
        {
            reads_to_save.insert(record.read_id);
        }
    }
    Ok(reads_to_save)
}

fn process_kraken_report_line(kraken_report: &str) -> Result<KrakenReportRecord> {
    let mut fields = kraken_report.split('\t');

    let percent_field = fields
        .next()
        .ok_or_else(|| eyre!("Missing percent field in the kraken report file"))?;
    let fragments_clade_rooted_field = fields
        .next()
        .ok_or_else(|| eyre!("Missing fragments clade rooted field in the kraken report file"))?;
    let fragments_taxon_field = fields
        .next()
        .ok_or_else(|| eyre!("Missing fragments taxon field in the kraken report file"))?;
    let rank_field = fields
        .next()
        .ok_or_else(|| eyre!("Missing rank field in the kraken report file"))?;
    let taxon_id_field = fields
        .next()
        .ok_or_else(|| eyre!("Missing taxon ID field in the kraken report file"))?;
    let name_field = fields
        .next()
        .ok_or_else(|| eyre!("Missing taxon name field in the kraken report file"))?;
    if fields.next().is_none() {
        let percent = percent_field
            .trim()
            .parse::<f32>()
            .wrap_err_with(|| format!("Error parsing percent value: '{}'", percent_field))?;

        let fragments_clade_rooted = fragments_clade_rooted_field
            .trim()
            .parse::<i32>()
            .wrap_err_with(|| {
                format!(
                    "Error parsing fragments clade rooted: '{}'",
                    fragments_clade_rooted_field
                )
            })?;

        let fragments_taxon = fragments_taxon_field
            .trim()
            .parse::<i32>()
            .wrap_err_with(|| {
                format!("Error parsing fragments taxon: '{}'", fragments_taxon_field)
            })?;

        let taxon_id = taxon_id_field
            .trim()
            .parse::<i32>()
            .wrap_err_with(|| format!("Error parsing taxon ID: '{}'", taxon_id_field))?;

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

pub fn build_tree_from_kraken_report(
    taxon_to_save: &[i32],
    report_path: &PathBuf,
) -> Result<(Vec<Tree>, HashMap<i32, usize>)> {
    debug!("Building taxonomic tree from kraken report");
    // will store the tree
    let mut nodes = Vec::new();
    // taxonid -> index in the nodes vector
    let mut taxon_map = HashMap::new();

    let report_file = fs::File::open(report_path)
        .wrap_err_with(|| format!("Failed to open kraken report file: {:?}", report_path))?;

    let reader = BufReader::new(report_file);
    let mut prev_index = None;

    for line in reader.lines() {
        let line = line.wrap_err("Error reading kraken report line")?;
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
            continue;
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
        if taxon_to_save.contains(&record.taxon_id) {
            taxon_map.insert(record.taxon_id, curr_index);
        }
    }

    let missing_taxon_ids = taxon_to_save
        .iter()
        .filter(|txid| !taxon_map.contains_key(txid))
        .collect::<Vec<&i32>>();

    if !missing_taxon_ids.is_empty() {
        bail!(
            "Failed to find taxon IDs in kraken report: {:?}",
            missing_taxon_ids
        );
    }

    Ok((nodes, taxon_map))
}

pub fn extract_parents(
    taxon_map: &HashMap<i32, usize>,
    nodes: &[Tree],
    taxon_id: i32,
) -> Result<Vec<i32>> {
    // Backtracking traversal from the given taxon_id to the root
    let start_index = taxon_map
        .get(&taxon_id)
        .ok_or_else(|| eyre!("Taxon ID {} not found in taxonomy map", taxon_id))?;

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

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write;
    use tempfile::tempdir;

    // kraken output tests
    #[test]
    fn test_process_kraken_output_line_valid() {
        let line = "C\tread_123\t1337\t150\t0:1 1:10";
        let result = process_kraken_output_line(line).unwrap();
        assert!(result.is_classified);
        assert_eq!(result.read_id, b"read_123");
        assert_eq!(result.taxon_id, 1337);
        assert_eq!(result.length, "150");
        assert_eq!(result.lca_map, "0:1 1:10");
    }

    #[test]
    fn test_process_kraken_output_line_unclassified() {
        let line = "U\tread_123\t1337\t150\t0:1 1:10";
        let result = process_kraken_output_line(line).unwrap();
        assert!(!result.is_classified);
        assert_eq!(result.read_id, b"read_123");
        assert_eq!(result.taxon_id, 1337);
        assert_eq!(result.length, "150");
        assert_eq!(result.lca_map, "0:1 1:10");
    }

    #[test]
    fn test_process_kraken_output_line_missing_fields() {
        let line = "";
        assert!(process_kraken_output_line(line).is_err());
        let line = "C";
        assert!(process_kraken_output_line(line).is_err());
        let line = "C\tread_123";
        assert!(process_kraken_output_line(line).is_err());
        let line = "C\tread_123\t1337";
        assert!(process_kraken_output_line(line).is_err());
        let line = "C\tread_123\t1337\t150";
        assert!(process_kraken_output_line(line).is_err());
    }

    #[test]
    fn test_process_kraken_output_line_invalid_taxon_id() {
        let line = "C\tread_123\tAAAA\t150\t0:1 1:10";
        assert!(process_kraken_output_line(line).is_err());
    }

    #[test]
    fn test_process_kraken_output_line_too_many_fields() {
        let line = "C\tread_123\t9606\t150\t0:1 1:10\tishouldntbehere";
        assert!(process_kraken_output_line(line).is_err());
    }

    #[test]
    fn test_process_kraken_output_include_mode() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("kraken_output.txt");
        let test_data = "\
        C\tread_1\t1337\t150\t0:1 1:10
        C\tread_2\t2\t150\t0:1 1:10
        C\tread_3\t1337\t150\t0:1 1:10
        U\tread_4\t0\t150\t0:1 1:10";
        let mut file = File::create(&file_path).unwrap();
        file.write_all(test_data.as_bytes()).unwrap();
        let taxon_ids_to_save = vec![1337];
        let result = process_kraken_output(&file_path, false, &taxon_ids_to_save).unwrap();
        assert_eq!(result.len(), 2);
        assert!(result.contains(b"read_1".as_slice()));
        assert!(result.contains(b"read_3".as_slice()));
        assert!(!result.contains(b"read_2".as_slice()));
        assert!(!result.contains(b"read_4".as_slice()));
    }

    #[test]
    fn test_process_kraken_output_exclude_mode() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("kraken_output.txt");
        let test_data = "\
        C\tread_1\t1337\t150\t0:1 1:10
        C\tread_2\t2\t150\t0:1 1:10
        C\tread_3\t1337\t150\t0:1 1:10
        U\tread_4\t0\t150\t0:1 1:10";
        let mut file = File::create(&file_path).unwrap();
        file.write_all(test_data.as_bytes()).unwrap();
        let taxon_ids_to_save = vec![1337];
        let result = process_kraken_output(&file_path, true, &taxon_ids_to_save).unwrap();
        assert_eq!(result.len(), 2);
        assert!(!result.contains(b"read_1".as_slice()));
        assert!(!result.contains(b"read_3".as_slice()));
        assert!(result.contains(b"read_2".as_slice()));
        assert!(result.contains(b"read_4".as_slice()));
    }

    #[test]
    fn test_process_kraken_output_empty_taxon_ids() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("kraken_output.txt");
        let test_data = "\
        C\tread_1\t1337\t150\t0:1 1:10
        C\tread_2\t2\t150\t0:1 1:10";
        let mut file = File::create(&file_path).unwrap();
        file.write_all(test_data.as_bytes()).unwrap();
        let result = process_kraken_output(&file_path, false, &[]).unwrap();
        assert_eq!(result.len(), 0);
        let result = process_kraken_output(&file_path, true, &[]).unwrap();
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_process_kraken_output_file_not_found() {
        let nonexistent_path = PathBuf::from("nonexistent_file.txt");
        let result = process_kraken_output(&nonexistent_path, false, &[1337]);
        assert!(result.is_err());
    }

    #[test]
    fn test_process_kraken_output_invalid_line() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("kraken_output.txt");
        let test_data = "\
        C\tread_1\t1337\t150\t0:1 1:10
        im_very_invalid
        C\tread_3\t1337\t150\t0:1 1:10";
        let mut file = File::create(&file_path).unwrap();
        file.write_all(test_data.as_bytes()).unwrap();
        let result = process_kraken_output(&file_path, false, &[1337]);
        assert!(result.is_err());
    }

    // kraken report tests
    #[test]
    fn test_process_kraken_report_line_valid() {
        let line = "10.77\t100\t50\tS\t1337\t  Homo sapiens";
        let result = process_kraken_report_line(line).unwrap();
        assert_eq!(result.percent, 10.77);
        assert_eq!(result.fragments_clade_rooted, 100);
        assert_eq!(result.fragments_taxon, 50);
        assert_eq!(result.rank, "S");
        assert_eq!(result.taxon_id, 1337);
        assert_eq!(result.level, 1);
        assert_eq!(result.name, "  Homo sapiens");
    }

    #[test]
    fn test_process_kraken_report_line_level() {
        let line = "5.2\t80\t30\tG\t1234\t      Escherichia";
        let result = process_kraken_report_line(line).unwrap();
        assert_eq!(result.level, 3);
        assert_eq!(result.name, "      Escherichia");
    }

    #[test]
    fn test_process_kraken_report_line_level_none() {
        let line = "90.0\t1000\t900\tD\t2\tBacteria";
        let result = process_kraken_report_line(line).unwrap();
        assert_eq!(result.level, 0);
        assert_eq!(result.name, "Bacteria");
    }

    #[test]
    fn test_process_kraken_report_line_missing_fields() {
        let line = "";
        assert!(process_kraken_report_line(line).is_err());
        let line = "10.5";
        assert!(process_kraken_report_line(line).is_err());
        let line = "10.5\t100";
        assert!(process_kraken_report_line(line).is_err());
        let line = "10.5\t100\t50";
        assert!(process_kraken_report_line(line).is_err());
        let line = "10.5\t100\t50\tS";
        assert!(process_kraken_report_line(line).is_err());
        let line = "10.5\t100\t50\tS\t1337";
        assert!(process_kraken_report_line(line).is_err());
    }

    #[test]
    fn test_process_kraken_report_line_invalid_fields() {
        let invalid_lines = vec![
            "not_a_number\t100\t50\tS\t1337\t  Homo sapiens",
            "10.5\tnot_a_number\t50\tS\t1337\t  Homo sapiens",
            "10.5\t100\tnot_a_number\tS\t1337\t  Homo sapiens",
            "10.5\t100\t50\tS\tnot_a_number\t  Homo sapiens",
        ];

        for line in invalid_lines {
            assert!(process_kraken_report_line(line).is_err());
        }
    }

    #[test]
    fn test_process_kraken_report_line_too_many_fields() {
        let line = "10.5\t100\t50\tS\t1337\t  Homo sapiens\textra_field";
        assert!(process_kraken_report_line(line).is_err());
    }

    // report parsing tests

    #[test]
    fn test_build_tree_from_kraken_report_valid() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("kraken_report.txt");
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
        let mut file = File::create(&file_path).unwrap();
        file.write_all(test_data.as_bytes()).unwrap();
        let taxon_to_save = vec![1386, 1239];
        let (nodes, taxon_map) = build_tree_from_kraken_report(&taxon_to_save, &file_path).unwrap();
        println!("{:?}", nodes);
        assert_eq!(nodes.len(), 10);

        // Check root
        assert_eq!(nodes[0].taxon_id, 1);
        assert_eq!(nodes[0].level_num, 0);
        assert_eq!(nodes[0].parent, None);
        assert_eq!(nodes[0].children, vec![1]);

        // Check Bacteria
        assert_eq!(nodes[2].taxon_id, 2);
        assert_eq!(nodes[2].level_num, 2);
        assert_eq!(nodes[2].parent, Some(1)); // Parent should be cellular organisms
        assert_eq!(nodes[2].children, vec![3]); // Children should be Terrabacteria group

        // Check Bacillota
        assert_eq!(nodes[4].taxon_id, 1239);
        assert_eq!(nodes[4].level_num, 4);
        assert_eq!(nodes[4].parent, Some(3)); // Parent should be Terrabacteria group
        assert_eq!(nodes[4].children, vec![5, 6]); // Children should be Bacilli and Bacilli2

        // Check that Bacilli2 and Bacilli are siblings (share same parent)
        assert_eq!(nodes[6].parent, Some(4));

        // Check taxon map
        assert_eq!(taxon_map.len(), 2);
        assert_eq!(taxon_map[&1386], 9);
        assert_eq!(taxon_map[&1239], 4);
    }

    #[test]
    fn test_build_tree_from_kraken_report_missing_taxon() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("kraken_report.txt");
        let test_data = "\
        21.36\t745591\t745591\tU\t0\tunclassified
        78.64\t2745487\t1646\tR\t1\troot
        78.58\t2743340\t1360\tR1\t131567\t  cellular organisms
        78.21\t2730479\t8458\tD\t2\t    Bacteria";
        let mut file = File::create(&file_path).unwrap();
        file.write_all(test_data.as_bytes()).unwrap();
        let taxon_to_save = vec![1386, 2];
        let result = build_tree_from_kraken_report(&taxon_to_save, &file_path);
        assert!(result.is_err());
    }

    #[test]
    fn test_build_tree_from_kraken_report_file_not_found() {
        let nonexistent_path = PathBuf::from("nonexistent_file.txt");
        let taxon_to_save = vec![1386];
        let result = build_tree_from_kraken_report(&taxon_to_save, &nonexistent_path);
        assert!(result.is_err());
    }

    #[test]
    fn test_build_tree_from_kraken_report_invalid_line() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("kraken_report.txt");
        let test_data = "\
        21.36\t745591\t745591\tU\t0\tunclassified
        78.64\t2745487\t1646\tR\t1\troot
        IM_AN_INVALID_LINE(((>?<???
        78.58\t2743340\t1360\tR1\t131567\t  cellular organisms";
        let mut file = File::create(&file_path).unwrap();
        file.write_all(test_data.as_bytes()).unwrap();
        let taxon_to_save = vec![131567];
        let result = build_tree_from_kraken_report(&taxon_to_save, &file_path);
        assert!(result.is_err());
    }

    #[test]
    fn test_extract_parents_valid() {
        let nodes = vec![
            Tree::new(1, 0, None),
            Tree::new(10, 1, Some(0)),
            Tree::new(20, 2, Some(1)),
            Tree::new(30, 3, Some(2)),
        ];
        let mut tree = nodes.clone();
        tree[0].children.push(1);
        tree[1].children.push(2);
        tree[2].children.push(3);
        let mut taxon_map = HashMap::new();
        taxon_map.insert(1, 0);
        taxon_map.insert(10, 1);
        taxon_map.insert(20, 2);
        taxon_map.insert(30, 3);
        let parents = extract_parents(&taxon_map, &tree, 30).unwrap();
        assert_eq!(parents, vec![30, 20, 10, 1]);
        let parents = extract_parents(&taxon_map, &tree, 20).unwrap();
        assert_eq!(parents, vec![20, 10, 1]);
        let parents = extract_parents(&taxon_map, &tree, 10).unwrap();
        assert_eq!(parents, vec![10, 1]);
        let parents = extract_parents(&taxon_map, &tree, 1).unwrap();
        assert_eq!(parents, vec![1]);
    }

    #[test]
    fn test_extract_parents_taxon_not_found() {
        let nodes = vec![Tree::new(1, 0, None)];
        let taxon_map = HashMap::new();
        let result = extract_parents(&taxon_map, &nodes, 99999);
        assert!(result.is_err());
    }

    #[test]
    fn test_extract_parents_invalid_parent_index() {
        let nodes = vec![Tree::new(1, 0, None), Tree::new(10, 1, Some(99999))];
        let mut taxon_map = HashMap::new();
        taxon_map.insert(10, 1);
        let result = extract_parents(&taxon_map, &nodes, 10);
        assert!(result.is_err());
    }

    #[test]
    fn test_extract_parents_empty_tree() {
        let nodes: Vec<Tree> = Vec::new();
        let taxon_map = HashMap::new();
        let result = extract_parents(&taxon_map, &nodes, 0);
        assert!(result.is_err());
    }

    #[test]
    fn test_extract_children_valid() {
        let mut nodes = vec![
            Tree::new(1, 0, None),
            Tree::new(10, 1, Some(0)),
            Tree::new(20, 1, Some(0)),
            Tree::new(30, 2, Some(1)),
            Tree::new(40, 2, Some(1)),
        ];
        nodes[0].children = vec![1, 2];
        nodes[1].children = vec![3, 4];
        let mut result = Vec::new();
        extract_children(&nodes, 0, &mut result).unwrap();
        assert_eq!(result, vec![30, 40, 10, 20, 1]);
        let mut result = Vec::new();
        extract_children(&nodes, 1, &mut result).unwrap();
        assert_eq!(result, vec![30, 40, 10]);
        let mut result = Vec::new();
        extract_children(&nodes, 3, &mut result).unwrap();
        assert_eq!(result, vec![30]);
    }

    #[test]
    fn test_extract_children_invalid_start_index() {
        let nodes = vec![Tree::new(1, 0, None)];
        let mut result = Vec::new();
        let out_of_bounds_idx = 999999;
        let extraction_result = extract_children(&nodes, out_of_bounds_idx, &mut result);
        assert!(extraction_result.is_err());
    }

    #[test]
    fn test_extract_children_invalid_child_index() {
        let mut nodes = vec![Tree::new(1, 0, None), Tree::new(10, 1, Some(0))];
        nodes[0].children.push(999);
        let mut result = Vec::new();
        let extraction_result = extract_children(&nodes, 0, &mut result);
        assert!(extraction_result.is_err());
    }

    #[test]
    fn test_extract_children_empty_tree() {
        let nodes: Vec<Tree> = Vec::new();
        let mut result = Vec::new();
        let extraction_result = extract_children(&nodes, 0, &mut result);
        assert!(extraction_result.is_err());
    }
}
