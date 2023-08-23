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
    /// Creates a new empty `FastaIndex`.
    pub fn new() -> Self {
        Self {
            entries: HashMap::new(),
        }
    }
    /// Inserts an `IndexEntry` into the `FastaIndex`.
    pub fn insert(&mut self, entry: IndexEntry) {
        self.entries.insert(entry.name.clone(), entry);
    }
    /// Creates a new `FastaIndex` from a `Read` object.
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
    /// Creates a new `FastaIndex` from a file path.
    pub fn from_filepath(path: &str) -> Result<Self> {
        let file = File::open(path)?;
        Self::from_reader(file)
    }
    /// Returns a reference to the `IndexEntry` corresponding to the given name.
    pub fn get(&self, name: &str) -> Option<&IndexEntry> {
        self.entries.get(name)
    }
    /// Returns a reference to the internal `HashMap` of entries.
    pub fn get_entries(&self) -> &HashMap<String, IndexEntry> {
        &self.entries
    }
}

#[cfg(test)]
mod testing {
    use crate::FastaIndex;
    use anyhow::Result;
    const TEST_FASTA_INDEX: &str = "example_data/example.fa.fai";

    #[test]
    fn build_index() -> Result<()> {
        let index = FastaIndex::from_filepath(TEST_FASTA_INDEX)?;
        assert_eq!(index.get_entries().len(), 2);
        Ok(())
    }
}
