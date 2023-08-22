use crate::IndexEntry;
use anyhow::Result;
use hashbrown::HashMap;
use std::{fs::File, io::Read};

/// A FASTA index.
///
/// This struct builds a map of FASTA entry names to their corresponding
/// `IndexEntry` structs.
#[derive(Debug)]
pub struct FastaIndex {
    entries: HashMap<String, IndexEntry>,
}
impl FastaIndex {
    pub fn new() -> Self {
        Self {
            entries: HashMap::new(),
        }
    }
    pub fn insert(&mut self, entry: IndexEntry) {
        self.entries.insert(entry.name.clone(), entry);
    }
    pub fn from_reader<R: Read>(reader: R) -> Result<Self> {
        let mut csv_reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(reader);
        let mut index = Self::new();
        for record in csv_reader.deserialize() {
            let record: IndexEntry = record?;
            index.insert(record);
        }
        Ok(index)
    }
    pub fn from_filepath(path: &str) -> Result<Self> {
        let file = File::open(path)?;
        Self::from_reader(file)
    }
    pub fn get(&self, name: &str) -> Option<&IndexEntry> {
        self.entries.get(name)
    }
}
