//! # faiquery
//!
//! `faiquery` is a library for querying FASTA files using the FAI index format.
//! It is designed to be fast and memory efficient, and is suitable for use in
//! high-throughput applications.
//!
//! It keeps a memory-mapped index of the FASTA file, and uses this to query
//! the file on demand using interval queries.
//!
//! ## Mutability
//!
//! The `IndexedFasta` has the option of keeping an internal buffer for
//! reading. This buffer is reused for all queries, and is cleared after each
//! query. This is the default behaviour.
//!
//! It will remove all newlines from the resulting queries, and return the
//! resulting sequence as a `&[u8]`.
//!
//! However, if you need to keep the newlines, or if you need to keep the
//! memory usage low then you can use the `query_buffer` method instead. This
//! will return a `&[u8]` directly from the memory map and avoid copying the
//! sequence into a buffer.
//! This will not remove newlines from the resulting sequence.
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
//! Let's show the default behavior which includes keeping an internal buffer
//! and requires a mutable `IndexedFasta` object.
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
//!
//! ### Querying the FASTA file immutably
//!
//! Let's now show the immutable behavior which does not keep an internal
//! buffer and does not require a mutable `IndexedFasta` object.
//!
//! ```rust
//! use faiquery::{FastaIndex, IndexedFasta};
//! use anyhow::Result;
//!
//! let index = FastaIndex::from_filepath("example_data/example.fa.fai")
//!     .expect("Could not read index file");
//! let faidx = IndexedFasta::new(index, "example_data/example.fa")
//!     .expect("Could not read FASTA file");
//!
//! // Query the first 10 bases of chr1
//! let seq = faidx.query_buffer("chr1", 0, 10).unwrap();
//! assert_eq!(seq, b"ACCTACGATC");
//!
//! // Query the first 10 bases of chr2
//! let seq = faidx.query_buffer("chr2", 0, 10).unwrap();
//! assert_eq!(seq, b"TTTTGATCGA");
//!
//! // Query the first 40 bases of chr1
//! let seq = faidx.query_buffer("chr1", 0, 40).unwrap();
//!
//! // The resulting sequence is 41 characters long
//! // This is because 1 newline is included
//! assert_eq!(seq.len(), 41);
//!
//! // The resulting sequence has 1 newline
//! let num_newlines = seq.iter().filter(|&&b| b == b'\n').count();
//! assert_eq!(num_newlines, 1);
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
    fn buffered_usage() -> Result<()> {
        let index = FastaIndex::from_filepath(TEST_FASTA_INDEX)?;
        let faidx = IndexedFasta::new(index, TEST_FASTA)?;
        let seq = faidx.query_buffer("chr1", 0, 10)?;
        assert_eq!(seq, b"ACCTACGATC");
        let seq = faidx.query_buffer("chr2", 0, 10)?;
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
    fn interval_over_newline_buffer() -> Result<()> {
        let index = FastaIndex::from_filepath(TEST_FASTA_INDEX)?;
        let faidx = IndexedFasta::new(index, TEST_FASTA)?;
        let seq = faidx.query_buffer("chr1", 20, 30)?;
        assert_eq!(seq, b"AGCTAGCT\nCA");
        let seq = faidx.query_buffer("chr2", 20, 30)?;
        assert_eq!(seq, b"CGCGCGGC\nCA");
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
    fn interval_overexted_left_buffered() -> Result<()> {
        let index = FastaIndex::from_filepath(TEST_FASTA_INDEX)?;
        let faidx = IndexedFasta::new(index, TEST_FASTA)?;
        let seq = faidx.query_buffer("chr1", 130, 150);
        assert!(seq.is_err());
        Ok(())
    }

    #[test]
    fn interval_overextend_right_buffered() -> Result<()> {
        let index = FastaIndex::from_filepath(TEST_FASTA_INDEX)?;
        let faidx = IndexedFasta::new(index, TEST_FASTA)?;
        let seq = faidx.query_buffer("chr1", 100, 150);
        assert!(seq.is_err());
        Ok(())
    }

    #[test]
    fn interval_overextend_start_eq_buffered() -> Result<()> {
        let index = FastaIndex::from_filepath(TEST_FASTA_INDEX)?;
        let faidx = IndexedFasta::new(index, TEST_FASTA)?;
        let seq = faidx.query_buffer("chr1", 112, 113);
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
    fn interval_overextend_left_unbounded_buffered() -> Result<()> {
        let index = FastaIndex::from_filepath(TEST_FASTA_INDEX)?;
        let faidx = IndexedFasta::new(index, TEST_FASTA)?;
        let seq = faidx.query_buffer_unbounded("chr1", 130, 150);
        assert!(seq.is_err());
        Ok(())
    }

    #[test]
    fn interval_overextend_right_unbounded_buffered() -> Result<()> {
        let index = FastaIndex::from_filepath(TEST_FASTA_INDEX)?;
        let faidx = IndexedFasta::new(index, TEST_FASTA)?;
        let seq = faidx.query_buffer_unbounded("chr1", 100, 150)?;
        assert_eq!(seq.len(), 13);
        assert_eq!(seq.iter().filter(|&&b| b == b'\n').count(), 1);
        Ok(())
    }

    #[test]
    fn interval_overextend_right_unbounded_start_eq_buffered() -> Result<()> {
        let index = FastaIndex::from_filepath(TEST_FASTA_INDEX)?;
        let faidx = IndexedFasta::new(index, TEST_FASTA)?;
        let seq = faidx.query_buffer_unbounded("chr1", 112, 150);
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
    fn missing_chr_buffered() -> Result<()> {
        let index = FastaIndex::from_filepath(TEST_FASTA_INDEX)?;
        let faidx = IndexedFasta::new(index, TEST_FASTA)?;
        let seq = faidx.query_buffer("chr3", 130, 150);
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
    fn malformed_interval_buffered() -> Result<()> {
        let index = FastaIndex::from_filepath(TEST_FASTA_INDEX)?;
        let faidx = IndexedFasta::new(index, TEST_FASTA)?;
        let seq = faidx.query_buffer("chr1", 130, 120);
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

    #[test]
    fn empty_interval_buffered() -> Result<()> {
        let index = FastaIndex::from_filepath(TEST_FASTA_INDEX)?;
        let faidx = IndexedFasta::new(index, TEST_FASTA)?;
        let seq = faidx.query_buffer("chr1", 130, 130);
        assert!(seq.is_err());
        Ok(())
    }
}
