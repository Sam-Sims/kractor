use lazy_static::lazy_static;
use std::sync::{Arc, Mutex};

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
    pub read_id: String,
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

lazy_static! {
    pub static ref TAXON_ID_COUNT: Arc<Mutex<usize>> = Arc::new(Mutex::new(0));
    pub static ref TAXON_IDS: Arc<Mutex<Vec<i32>>> = Arc::new(Mutex::new(Vec::new()));
    pub static ref TOTAL_READS: Arc<Mutex<usize>> = Arc::new(Mutex::new(0));
    pub static ref READS_TO_EXTRACT: Arc<Mutex<usize>> = Arc::new(Mutex::new(0));
}
