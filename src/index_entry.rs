use serde::{Deserialize, Serialize};

/// A FASTA index entry.
///
/// This struct represents a single entry in a FASTA index.
/// It contains the name of the entry, the length of the entry,
/// the offset of the entry in the FASTA file, and the line
/// width and line bases of the entry.
#[derive(Serialize, Deserialize, Debug)]
pub struct IndexEntry {
    pub name: String,
    pub length: usize,
    pub offset: usize,
    pub line_bases: usize,
    pub line_width: usize,
}
