use crate::{FastaIndex, IndexEntry};
use anyhow::{bail, Result};
use memmap2::Mmap;
use std::fs::File;

/// An indexed FASTA file.
///
/// This struct is used to query a FASTA file by name and position.
/// It uses a memory-mapped file to avoid loading the entire file into memory.
/// It requires a `FastaIndex` to be created first.
///
/// # Examples
///
/// ```
/// use faiquery::{FastaIndex, IndexedFasta};
///
/// let index = FastaIndex::from_filepath("example_data/example.fa.fai")
///     .expect("Could not read index file");
/// let mut faidx = IndexedFasta::new(index, "example_data/example.fa")
///     .expect("Could not read FASTA file");
///
/// // Query the first 10 bases of chr1
/// let seq = faidx.query("chr1", 0, 10).unwrap();
/// assert_eq!(seq, b"ACCTACGATC");
/// ```
#[derive(Debug)]
pub struct IndexedFasta {
    index: FastaIndex,
    map: Mmap,
    buffer: Vec<u8>,
}
impl IndexedFasta {
    /// Create a new `IndexedFasta` from a `FastaIndex` and a file path.
    pub fn new(index: FastaIndex, path: &str) -> Result<Self> {
        let file = File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };
        let buffer = Vec::new();
        Ok(Self {
            index,
            map: mmap,
            buffer,
        })
    }

    /// Validate the start and end positions of a query interval.
    fn validate_interval(
        &self,
        entry: &IndexEntry,
        start: usize,
        end: usize,
        bounded: bool,
    ) -> Result<()> {
        if start > end {
            bail!("Start position must be less than end position");
        } else if start == end {
            bail!("Start and end positions must not be equal");
        } else if start >= entry.length {
            bail!("Start position must be less than sequence length");
        } else if bounded && end > entry.length {
            bail!("End position must be less than sequence length");
        }
        Ok(())
    }

    /// Query the FASTA file by name and position.
    ///
    /// The sequence is returned as a `&[u8]` slice but is not guaranteed to be valid UTF-8.
    /// It also removes all newline characters from the sequence slice.
    ///
    /// To query intervals that are potentially out of bounds, use `query_unbounded` and
    /// return a truncated sequence use `query_unbounded`.
    ///
    /// # Errors
    ///
    /// - Error if the query `name`is not found in the index.
    /// - Error if the `start` position is greater than the `end` position.
    /// - Error if the `start` position is equal to the `end` position.
    /// - Error if the `end` position is greater than the index sequence length.
    ///
    /// # Example
    ///
    /// ```
    /// use faiquery::{FastaIndex, IndexedFasta};
    ///
    /// let index = FastaIndex::from_filepath("example_data/example.fa.fai")
    ///    .expect("Could not read index file");
    /// let mut faidx = IndexedFasta::new(index, "example_data/example.fa")
    ///     .expect("Could not read FASTA file");
    ///
    /// // Query a region from the middle of chr1
    /// let seq = faidx.query("chr1", 50, 80).unwrap();
    /// assert_eq!(seq.len(), 30);
    ///
    /// // Overextend the query into chr1 (which is 112 bases long)
    /// let seq = faidx.query("chr1", 100, 120);
    /// assert!(seq.is_err());
    /// ```
    pub fn query(&mut self, name: &str, start: usize, end: usize) -> Result<&[u8]> {
        let entry = match self.index.get(name) {
            Some(entry) => entry,
            None => bail!("No entry found for {}", name),
        };
        self.validate_interval(entry, start, end, true)?;
        self.buffer.clear();
        let query_pos = QueryPosition::new(start, end, entry);
        let seq_slice = &self.map[query_pos.pos..query_pos.pos + query_pos.buffer_size];
        self.buffer.extend_from_slice(seq_slice);
        self.buffer.retain(|&c| c != b'\n');
        Ok(&self.buffer)
    }

    /// Query the FASTA file by name and position but do not copy to internal buffer.
    ///
    /// This will **not** remove newline characters from the sequence slice.
    ///
    /// This method is useful for memory-efficient operations where the sequence is not
    /// necessarily needed as a `&str` slice (since newlines are potentially included).
    ///
    /// # Errors
    ///
    /// - Error if the query `name` is not found in the index.
    /// - Error if the `start` position is greater than the `end` position.
    /// - Error if the `start` position is equal to the `end` position.
    /// - Error if the `end` position is greater than the index sequence length.
    ///
    /// ```
    /// use faiquery::{FastaIndex, IndexedFasta};
    ///
    /// let index = FastaIndex::from_filepath("example_data/example.fa.fai")
    ///    .expect("Could not read index file");
    /// let mut faidx = IndexedFasta::new(index, "example_data/example.fa")
    ///     .expect("Could not read FASTA file");
    ///
    /// // Query a region from the middle of chr1 (note the newline character is included)
    /// let seq = faidx.query_buffer("chr1", 50, 80).unwrap();
    /// assert_eq!(seq.len(), 31);
    /// assert!(seq.contains(&b'\n'));
    /// assert_eq!(seq.iter().filter(|&&c| c != b'\n').count(), 30);
    ///
    /// // Overextend the query into chr1 (which is 112 bases long)
    /// let seq = faidx.query_buffer("chr1", 100, 120);
    /// assert!(seq.is_err());
    /// ```
    pub fn query_buffer(&self, name: &str, start: usize, end: usize) -> Result<&[u8]> {
        let entry = match self.index.get(name) {
            Some(entry) => entry,
            None => bail!("No entry found for {}", name),
        };
        self.validate_interval(entry, start, end, true)?;
        let query_pos = QueryPosition::new(start, end, entry);
        let seq_slice = &self.map[query_pos.pos..query_pos.pos + query_pos.buffer_size];
        Ok(seq_slice)
    }

    /// Query the FASTA file by name and position.
    ///
    /// The sequence is returned as a `&[u8]` slice but is not guaranteed to be valid UTF-8.
    /// It also removes all newline characters from the sequence slice.
    ///
    /// This method will truncate the sequence if the `end` position is greater than the sequence length
    /// to avoid an error and only return the sequence up to the sequence length.
    ///
    /// # Errors
    ///
    /// - Error if the query `name`is not found in the index.
    /// - Error if the `start` position is greater than the `end` position.
    /// - Error if the `start` position is equal to the `end` position.
    ///
    /// # Example
    ///
    /// ```
    /// use faiquery::{FastaIndex, IndexedFasta};
    ///
    /// let index = FastaIndex::from_filepath("example_data/example.fa.fai")
    ///    .expect("Could not read index file");
    /// let mut faidx = IndexedFasta::new(index, "example_data/example.fa")
    ///     .expect("Could not read FASTA file");
    ///
    /// // Overextend the query into chr1 (which is 112 bases long)
    /// let seq = faidx.query("chr1", 100, 120);
    /// assert!(seq.is_err());
    ///
    /// // Overextend the query into chr1 but truncate the sequence
    /// // with `query_unbounded`
    /// let seq = faidx.query_unbounded("chr1", 100, 120).unwrap();
    /// assert_eq!(seq.len(), 12);
    /// ```
    pub fn query_unbounded(&mut self, name: &str, start: usize, end: usize) -> Result<&[u8]> {
        let entry = match self.index.get(name) {
            Some(entry) => entry,
            None => bail!("No entry found for {}", name),
        };
        self.validate_interval(entry, start, end, false)?;
        let end = if end > entry.length {
            entry.length
        } else {
            end
        };
        self.buffer.clear();
        let query_pos = QueryPosition::new(start, end, entry);
        let seq_slice = &self.map[query_pos.pos..query_pos.pos + query_pos.buffer_size];
        self.buffer.extend_from_slice(seq_slice);
        self.buffer.retain(|&c| c != b'\n');
        Ok(&self.buffer)
    }

    /// Query the FASTA file by name and position.
    ///
    /// The sequence is returned as a `&[u8]` slice but is not guaranteed to be valid UTF-8.
    /// This will **not** remove newline characters from the sequence slice.
    ///
    /// This method will truncate the sequence if the `end` position is greater than the sequence length
    /// to avoid an error and only return the sequence up to the sequence length.
    ///
    /// # Errors
    ///
    /// - Error if the query `name`is not found in the index.
    /// - Error if the `start` position is greater than the `end` position.
    /// - Error if the `start` position is equal to the `end` position.
    ///
    /// # Example
    ///
    /// ```
    /// use faiquery::{FastaIndex, IndexedFasta};
    ///
    /// let index = FastaIndex::from_filepath("example_data/example.fa.fai")
    ///    .expect("Could not read index file");
    /// let mut faidx = IndexedFasta::new(index, "example_data/example.fa")
    ///     .expect("Could not read FASTA file");
    ///
    /// // Overextend the query into chr1 (which is 112 bases long)
    /// let seq = faidx.query_buffer("chr1", 100, 120);
    /// assert!(seq.is_err());
    ///
    /// // Overextend the query into chr1 but truncate the sequence
    /// // with `query_unbounded`
    /// let seq = faidx.query_buffer_unbounded("chr1", 100, 120).unwrap();
    ///
    /// // The sequence is truncated to 13 characters
    /// assert_eq!(seq.len(), 13);
    ///
    /// // The sequence contains newline characters
    /// assert!(seq.contains(&b'\n'));
    ///
    /// // The sequence contains 12 non-newline characters
    /// assert_eq!(seq.iter().filter(|&&c| c != b'\n').count(), 12);
    /// ```
    pub fn query_buffer_unbounded(&self, name: &str, start: usize, end: usize) -> Result<&[u8]> {
        let entry = match self.index.get(name) {
            Some(entry) => entry,
            None => bail!("No entry found for {}", name),
        };
        self.validate_interval(entry, start, end, false)?;
        let end = if end > entry.length {
            entry.length
        } else {
            end
        };
        let query_pos = QueryPosition::new(start, end, entry);
        let seq_slice = &self.map[query_pos.pos..query_pos.pos + query_pos.buffer_size];
        Ok(seq_slice)
    }
}

/// A query position.
///
/// This struct is used to calculate the position of a query in a FASTA file.
/// It is used to calculate the offset and size of the query in the memory-mapped file.
struct QueryPosition {
    pub buffer_size: usize,
    pub pos: usize,
}
impl QueryPosition {
    pub fn new(start: usize, end: usize, entry: &IndexEntry) -> Self {
        let size = end - start;
        let row_pos = (start / entry.line_bases) * entry.line_width;
        let col_pos = start % entry.line_bases;
        let num_lines = (size + col_pos) / entry.line_bases;
        let buffer_size = size + num_lines;
        let pos = entry.offset + row_pos + col_pos;
        Self { buffer_size, pos }
    }
}
