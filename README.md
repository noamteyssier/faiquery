# faiquery

[![MIT licensed](https://img.shields.io/badge/license-MIT-blue.svg)](./LICENSE.md)
[![docs.rs](https://img.shields.io/docsrs/faiquery?color=green&label=docs.rs)](https://docs.rs/faiquery/latest/faiquery/)

perform interval queries on an indexed fasta file

## Description

This is a simple utility library to request interval queries
on an indexed fasta file.

Ths index is assumed [`samtools faidx`](http://www.htslib.org/doc/samtools-faidx.html).

The fasta file is memory-mapped and a single mutable buffer is kept
for the indexed fasta reader.
This buffer is used to return a slice reference to the query sequence
but allows for a non-contiguous sequence (i.e. can return sequences without newlines)

Note that because a mutable buffer is kept, this is not the best approach for
concurrent operations.

However, for single-threaded applications, this performs very well with low memory
overhead and runtime.

## Usage

```rust
use anyhow::Result;
use faiquery::{FastaIndex, IndexedFasta};

fn main() -> Result<()> {
    let index = FastaIndex::from_filepath("example_data/example.fa.fai")?;
    let mut faidx = IndexedFasta::new(index, "example_data/example.fa")?;

    // Query the first 10 bases of chr1
    let seq = faidx.query("chr1", 0, 10)?;
    assert_eq!(seq, b"ACCTACGATC");

    // Query the first 10 bases of chr2
    let seq = faidx.query("chr2", 0, 10)?;
    assert_eq!(seq, b"TTTTGATCGA");

    Ok(())
}
```

## Similar Approaches

There are other index fasta readers that are available,
here is a nonexhaustive list of those:

- [faimm](https://crates.io/crates/faimm)
- [noodles](https://crates.io/crates/noodles)
- [rust-bio](https://crates.io/crates/bio)

If you are looking for a powerful concurrent reader I recommend
using the [faimm](https://crates.io/crates/faimm) crate.
