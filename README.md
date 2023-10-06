# Butcher

Butcher is a read trimmer for single-end long-read sequencing (Nanopore) or paired-end Illumina sequencing. Much like Nanofilt, Porechop, Trimmomatic, 
and a host of other tools, Butcher can trim low-quality regions off the end of the reads, trim by adapter sequences, or remove poly-A and poly-G tracks
from sequencing reads. One nice feature is you can preview the changes in the terminal before running, which can be helpful to see how things would look
before committing to run. 

## Getting butcher

Downloads are available on the release page for Linux systems, or you can build it (see below).

## Documentation

```
A constrained use-case fastq trimmer

Usage: butcher [OPTIONS]

Options:
  -f, --fastq1 <FASTQ1>
          the first fastq file -- required
  -f, --fastq2 <FASTQ2>
          a second fastq file, for Illumina paired-end reads
  -o, --out-fastq1 <OUT_FASTQ1>
          output fastq file 1
  -o, --out-fastq2 <OUT_FASTQ2>
          output fastq file 2 -- if using paired-end reads
  -m, --minimum-remaining-read-size <MINIMUM_REMAINING_READ_SIZE>
          minimum remaining read size after trimming is complete -- reads shorter than this will be discarded [default: 10]
  -w, --window-min-qual-score <WINDOW_MIN_QUAL_SCORE>
          the minimum average quality score a window of nucleotides must have [default: 10]
  -w, --window-size <WINDOW_SIZE>
          trimming window size [default: 10]
  -t, --trim-poly-a
          enable poly-A tail trimming (seen in RNA-seq data)
  -t, --trim-poly-g
          trim a read after a poly-G tail is found (seen when sequencing off the end of an Illumina read with 2-color chemistry)
  -t, --trim-poly-x-length <TRIM_POLY_X_LENGTH>
          the length of the poly-X tail to trim (use with the poly-A or poly-G trimming) [default: 10]
  -t, --trim-poly-x-proportion <TRIM_POLY_X_PROPORTION>
          the proportion of bases that must be X to trim the read end [default: 0.9]
  -p, --primers <PRIMERS>
          primers to detect and remove (we'll make their reverse complement too), separated by commas
  -p, --primers-max-mismatch-distance <PRIMERS_MAX_MISMATCH_DISTANCE>
          the maximum mismatches allowed for primer trimming -- it's best if this is 1 or 2 [default: 1]
  -p, --primers-end-proportion <PRIMERS_END_PROPORTION>
          what proportion of the read ends can a primer be found in (front or back) -- if it's interior to this margin we drop the read(s) [default: 0.2]
  -p, --preview
          just display the reads and what we'd cut, don't actually write any output to disk
  -h, --help
          Print help

```

## Compiling

Rust nightly is required, mostly to avoid the crash text when you pipe (```|```) __butcher__ on the command line. Building with nightly is simple:
```
cargo build --release
```
From the command line. A __butcher__ artifact will be in the target/release/ folder.

