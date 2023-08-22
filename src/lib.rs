//! # faiquery
//!
//! `faiquery` is a library for querying FASTA files using the FAI index format.
//! It is designed to be fast and memory efficient, and is suitable for use in
//! high-throughput applications.
//!
//! It keeps a memory-mapped index of the FASTA file, and uses this to query
//! the file on demand using interval queries.
//! It keeps a single internal buffer for reading, and reuses this buffer for
//! all queries.
//! It removes all newlines from the resulting queries, and returns the
//! resulting sequence as a `&[u8]`.
//!
//! ## Example
//!
//! Here is an example fasta file:
//!
//! ### example.fa
//!
//! ```text
//! >chr1
//! ACCTACGATCGACTGATCGTAGCTAGCT
//! CATCGATCGTACGGACGATCGATCGGTT
//! CACACCGGGCATGACTGATCGGGGGCCC
//! ACGTGTGTGCAGCGCGCGGCGCGCGCGG
//! >chr2
//! TTTTGATCGATCGGCGGGCGCGCGCGGC
//! CAGATTCGGGCGCGATTATATATTAGCT
//! CGACGGCGACTCGAGCTACACGTCGGGC
//! GCGAGCGGGACGCGCGGCGCGCGCGGCC
//! AAAAAAATTTTTATATATTATTACGCGC
//! CGACTCAGTCGACTGGGGGCGCGCGCGC
//! AAACCACA
//! ```
//!
//! and its corresponding index file:
//!
//! ### example.fa.fai
//!
//! ```text
//! chr1	112	6	28	29
//! chr2	176	128	28	29
//! ```
//!
//! ### Querying the FASTA file
//!
//! ```rust
//! use faiquery::{FastaIndex, IndexedFasta};
//! use anyhow::Result;
//!
//! let index = FastaIndex::from_filepath("example_data/example.fa.fai")
//!     .expect("Could not read index file");
//! let mut faidx = IndexedFasta::new(index, "example_data/example.fa")
//!     .expect("Could not read FASTA file");
//!
//! // Query the first 10 bases of chr1
//! let seq = faidx.query("chr1", 0, 10).unwrap();
//! assert_eq!(seq, b"ACCTACGATC");
//!
//! // Query the first 10 bases of chr2
//! let seq = faidx.query("chr2", 0, 10).unwrap();
//! assert_eq!(seq, b"TTTTGATCGA");
//!
//! // Query the first 40 bases of chr1
//! let seq = faidx.query("chr1", 0, 40).unwrap();
//!
//! // The resulting sequence is 40 bases long
//! assert_eq!(seq.len(), 40);
//!
//! // The resulting sequence has no newlines
//! let num_newlines = seq.iter().filter(|&&b| b == b'\n').count();
//! assert_eq!(num_newlines, 0);
//! ```

mod fasta_index;
mod index_entry;
mod indexed_fasta;

/// The `FastaIndex` struct represents a FAI index file.
pub use fasta_index::FastaIndex;

/// The `IndexEntry` struct represents a single entry in a FAI index file.
pub use index_entry::IndexEntry;

/// The `IndexedFasta` struct represents a FASTA file that has been indexed
/// using the FAI format.
pub use indexed_fasta::IndexedFasta;

#[cfg(test)]
mod testing {
    use crate::{FastaIndex, IndexedFasta};
    use anyhow::Result;

    const TEST_FASTA: &str = "example_data/example.fa";
    const TEST_FASTA_INDEX: &str = "example_data/example.fa.fai";

    #[test]
    fn standard_usage() -> Result<()> {
        let index = FastaIndex::from_filepath(TEST_FASTA_INDEX)?;
        let mut faidx = IndexedFasta::new(index, TEST_FASTA)?;
        let seq = faidx.query("chr1", 0, 10)?;
        assert_eq!(seq, b"ACCTACGATC");
        let seq = faidx.query("chr2", 0, 10)?;
        assert_eq!(seq, b"TTTTGATCGA");
        Ok(())
    }

    #[test]
    fn interval_over_newline() -> Result<()> {
        let index = FastaIndex::from_filepath(TEST_FASTA_INDEX)?;
        let mut faidx = IndexedFasta::new(index, TEST_FASTA)?;
        let seq = faidx.query("chr1", 20, 30)?;
        assert_eq!(seq, b"AGCTAGCTCA");
        let seq = faidx.query("chr2", 20, 30)?;
        assert_eq!(seq, b"CGCGCGGCCA");
        Ok(())
    }

    #[test]
    fn interval_overextend_left() -> Result<()> {
        let index = FastaIndex::from_filepath(TEST_FASTA_INDEX)?;
        let mut faidx = IndexedFasta::new(index, TEST_FASTA)?;
        let seq = faidx.query("chr1", 130, 150);
        assert!(seq.is_err());
        Ok(())
    }

    #[test]
    fn interval_overextend_right() -> Result<()> {
        let index = FastaIndex::from_filepath(TEST_FASTA_INDEX)?;
        let mut faidx = IndexedFasta::new(index, TEST_FASTA)?;
        let seq = faidx.query("chr1", 100, 150);
        assert!(seq.is_err());
        Ok(())
    }

    #[test]
    fn interval_overextend_start_eq() -> Result<()> {
        let index = FastaIndex::from_filepath(TEST_FASTA_INDEX)?;
        let mut faidx = IndexedFasta::new(index, TEST_FASTA)?;
        let seq = faidx.query("chr1", 112, 113);
        assert!(seq.is_err());
        Ok(())
    }

    #[test]
    fn interval_overextend_left_unbounded() -> Result<()> {
        let index = FastaIndex::from_filepath(TEST_FASTA_INDEX)?;
        let mut faidx = IndexedFasta::new(index, TEST_FASTA)?;
        let seq = faidx.query_unbounded("chr1", 130, 150);
        assert!(seq.is_err());
        Ok(())
    }

    #[test]
    fn interval_overextend_right_unbounded() -> Result<()> {
        let index = FastaIndex::from_filepath(TEST_FASTA_INDEX)?;
        let mut faidx = IndexedFasta::new(index, TEST_FASTA)?;
        let seq = faidx.query_unbounded("chr1", 100, 150)?;
        assert_eq!(seq.len(), 12);
        Ok(())
    }

    #[test]
    fn interval_overextend_right_unbounded_start_eq() -> Result<()> {
        let index = FastaIndex::from_filepath(TEST_FASTA_INDEX)?;
        let mut faidx = IndexedFasta::new(index, TEST_FASTA)?;
        let seq = faidx.query_unbounded("chr1", 112, 150);
        assert!(seq.is_err());
        Ok(())
    }

    #[test]
    fn missing_chr() -> Result<()> {
        let index = FastaIndex::from_filepath(TEST_FASTA_INDEX)?;
        let mut faidx = IndexedFasta::new(index, TEST_FASTA)?;
        let seq = faidx.query("chr3", 130, 150);
        assert!(seq.is_err());
        Ok(())
    }

    #[test]
    fn malformed_interval() -> Result<()> {
        let index = FastaIndex::from_filepath(TEST_FASTA_INDEX)?;
        let mut faidx = IndexedFasta::new(index, TEST_FASTA)?;
        let seq = faidx.query("chr1", 130, 120);
        assert!(seq.is_err());
        Ok(())
    }

    #[test]
    fn empty_interval() -> Result<()> {
        let index = FastaIndex::from_filepath(TEST_FASTA_INDEX)?;
        let mut faidx = IndexedFasta::new(index, TEST_FASTA)?;
        let seq = faidx.query("chr1", 130, 130);
        assert!(seq.is_err());
        Ok(())
    }
}
