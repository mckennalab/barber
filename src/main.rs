#![feature(unix_sigpipe)]

extern crate colored;
extern crate clap;
extern crate flate2;
extern crate core;

mod trimmers;
mod sift4;
mod primers;

use std::cmp::{max, min};
use std::io;
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use colored::Colorize;
use clap::Parser;
use flate2::Compression;
use flate2::write::GzEncoder;
use crate::trimmers::{BackTrimmer, FastqTrimmer, FrontBackTrimmer, PolyXTrimmer, PrimerTrimmer};
use log::{warn};

/// Simple program to greet a person
#[derive(Parser, Debug)]
struct Args {
    /// input fastq file 1
    #[arg(short, long)]
    fastq1: Option<String>,

    /// input fastq file 2
    #[arg(short, long)]
    fastq2: Option<String>,

    /// output fastq file 1
    #[arg(short, long)]
    out_fastq1: Option<String>,

    /// output fastq file 2
    #[arg(short, long)]
    out_fastq2: Option<String>,

    /// minimum remaining read size
    #[arg(short, long, default_value_t = 10)]
    minimum_remaining_read_size: usize,

    /// the minimum average score a window of nucleotides must have
    #[arg(short, long, default_value_t = 10)]
    window_min_qual_score: u8,

    /// trimming window size
    #[arg(short, long, default_value_t = 10)]
    window_size: u8,

    /// trim a read after a poly-A tail is foun (RNA-seq)
    #[arg(short, long, default_value_t = false)]
    trim_poly_a: bool,

    /// trim a read after a poly-G tail is found (sequenced off the end of an Illumina read)
    #[arg(short, long, default_value_t = false)]
    trim_poly_g: bool,

    /// the length of the poly-X tral to trim (above)
    #[arg(short, long, default_value_t = 10)]
    trim_poly_x_length: usize,

    /// the length of the poly-X tral to trim (above)
    #[arg(short, long, default_value_t = 0.9)]
    trim_poly_x_proportion: f64,

    /// primers to detect and remove, separated by commas
    #[arg(short, long)]
    primers: Option<String>,

    /// primers max mismatch distance
    #[arg(short, long, default_value_t = 1)]
    primers_max_mismatch_distance: u8,

    /// what proportion of the read ends can a primer be found in -- if it's interior to this margin we drop the read(s)
    #[arg(short, long, default_value_t = 0.2)]
    primers_end_proportion: f64,

    /// just display the reads and what we'd cut
    #[arg(short, long, default_value_t = false)]
    preview: bool,
}

/// a simple FASTQ record with name, sequence, and quality
pub struct FastqRecord {
    pub name: Vec<u8>,
    pub seq: Vec<u8>,
    pub quals: Vec<u8>,
}

/// an input decoder for our gzipped FASTQ file
struct FastqInputFile {
    decoder: BufReader<MultiGzDecoder<File>>,
}

impl FastqInputFile {
    pub fn new(path: &str) -> Result<FastqInputFile, io::Error> {
        let file = File::open(path)?;
        let decoder = io::BufReader::new(MultiGzDecoder::new(file));
        Ok(FastqInputFile { decoder })
    }
}

impl Iterator for FastqInputFile {
    type Item = FastqRecord;

    fn next(&mut self) -> Option<FastqRecord> {
        let mut name = String::new();
        match self.decoder.read_line(&mut name) {
            Ok(_) => {
            }
            Err(_e) => {
                warn!("Error reading sequence line for unnamed read");
                return None
            },
        }
        let mut name = name.into_bytes();
        if name.len() == 0 {
            return None;

        }
        assert_eq!(name[0], b'@');
        name.pop(); // drop endline

        let mut seq = String::new();
        match self.decoder.read_line(&mut seq) {
            Ok(_) => {}
            Err(_e) => {
                warn!("Error reading sequence line for read {}", String::from_utf8(name).unwrap());
                return None
            },
        }
        seq.pop(); // drop endline
        let mut _orient = String::new();
        match self.decoder.read_line(&mut _orient) {
            Ok(_) => {}
            Err(_e) => {
                warn!("Error reading orientation line for read {}", String::from_utf8(name).unwrap());
                return None
            },
        }
        let mut quals = String::new();
        match self.decoder.read_line(&mut quals) {
            Ok(_) => {}
            Err(_e) => {
                warn!("Error reading quals line for read {}", String::from_utf8(name).unwrap());
                return None
            },
        }
        quals.pop(); // drop endline

        Some(FastqRecord { name, seq: seq.into_bytes(), quals: quals.into_bytes() })
    }
}

#[unix_sigpipe = "sig_dfl"]
fn main() {
    simple_logger::init_with_level(log::Level::Warn).unwrap();

    let args = Args::parse();

    assert!(args.preview ^ args.out_fastq1.is_some(), "Either preview mode or output files need to be set");

    let mut out_fastq1 = setup_compressed_file(&args.out_fastq1);
    let mut out_fastq2 = setup_compressed_file(&args.out_fastq2);

    let mut reader = FastqInputFile::new(&args.fastq1.unwrap()).expect("invalid path/file for fastq1");

    let mut cutters: Vec<Box<dyn FastqTrimmer>> = Vec::new();

    if args.primers.is_some() {
        let primers = args.primers.unwrap();
        let primers: Vec<Vec<u8>> = primers.split(",").map(|i|i.as_bytes().to_vec()).collect();
        println!("Using primers: {}", &primers.clone().into_iter().map(|i|String::from_utf8(i).unwrap()).collect::<Vec<String>>().join(", "));

        cutters.push(Box::new(PrimerTrimmer::new(&primers, &args.primers_max_mismatch_distance, &args.primers_end_proportion)));
    }

    if args.trim_poly_a {
        cutters.push(Box::new(PolyXTrimmer {
            window_size: args.trim_poly_x_length.clone(),
            minimum_g_proportion: args.trim_poly_x_proportion.clone(),
            bases: vec![b'A', b'a'],
        }));
    }
    if args.trim_poly_g {
        cutters.push(Box::new(PolyXTrimmer {
            window_size: args.trim_poly_x_length.clone(),
            minimum_g_proportion: args.trim_poly_x_proportion.clone(),
            bases: vec![b'G', b'g'],
        }));
    }

    if args.fastq2.is_some() {
        cutters.push(Box::new(BackTrimmer { window_size: args.window_size.clone(), window_min_qual_score: args.window_min_qual_score, qual_score_base: 32 }));
        let mut reader2 = FastqInputFile::new(&args.fastq2.unwrap()).expect("invalid path/file for fastq2");
        paired_end(&mut reader, &mut reader2, &mut out_fastq1, &mut out_fastq2, &cutters, &args.minimum_remaining_read_size, &args.preview);
    } else {
        cutters.push(Box::new(FrontBackTrimmer { window_size: args.window_size.clone(), window_min_qual_score: args.window_min_qual_score, qual_score_base: 32 }));
        single_end(&mut reader, &mut out_fastq1, &cutters, &args.minimum_remaining_read_size, &args.preview);
    }
}

fn max_min_pair(pair_one: &(bool, usize, usize), pair_two: &(bool, usize, usize)) -> (bool, usize, usize) {
    if !pair_one.0 || !pair_two.0 {
        (false, 0, 0)
    } else {
        (true, max(pair_one.1, pair_two.1), min(pair_one.2, pair_two.2))
    }
}

fn single_end(reader1: &mut FastqInputFile,
              out_fastq: &mut BufWriter<GzEncoder<Box<dyn Write>>>,
              cutters: &Vec<Box<dyn FastqTrimmer>>,
              minimum_remaining_read_size: &usize,
              preview: &bool) {

    while let Some(read1) = reader1.next() {
        let mut base_cuts = (true, 0, read1.seq.len());
        for cutter in cutters {
            let cut = cutter.trim(&read1);
            base_cuts = max_min_pair(&base_cuts, &cut);
        }

        if base_cuts.0 && &base_cuts.2 - &base_cuts.1 > *minimum_remaining_read_size {
            match *preview {
                true => {
                    print_read(read1.name.as_slice(), read1.seq.as_slice(), read1.quals.as_slice(), &base_cuts.1, &base_cuts.2)
                }
                false => {
                    write_read(out_fastq, &slice_read(&read1, &base_cuts.1, &base_cuts.2)).expect("Unable to write to output file 1.");
                }
            }
        }
    }
    out_fastq.flush().expect("Unable to flush output fastq file.");
}

fn paired_end(reader1: &mut FastqInputFile,
              reader2: &mut FastqInputFile,
              out_fastq1: &mut BufWriter<GzEncoder<Box<dyn Write>>>,
              out_fastq2: &mut BufWriter<GzEncoder<Box<dyn Write>>>,
              cutters: &Vec<Box<dyn FastqTrimmer>>,
              minimum_remaining_read_size: &usize,
              preview: &bool) {

    while let Some(read1) = reader1.next() {
        let read2 = match reader2.next() {
            None => {panic!("Reads in fastq1 and fastq2 are not paired, at read1 {}",String::from_utf8(read1.name).unwrap())}
            Some(x) => {x}
        };

        let mut base_cuts_read1 = (true, 0, read1.seq.len());
        let mut base_cuts_read2 = (true, 0, read2.seq.len());

        for cutter in cutters {
            let cut = cutter.trim(&read1);
            base_cuts_read1 = max_min_pair(&base_cuts_read1, &cut);

            let cut = cutter.trim(&read2);
            base_cuts_read2 = max_min_pair(&base_cuts_read2, &cut);
        }

        if ((base_cuts_read1.0 && base_cuts_read2.0) &&
            &base_cuts_read1.2 - &base_cuts_read1.1 > *minimum_remaining_read_size) &&
            (&base_cuts_read2.2 - &base_cuts_read2.1 > *minimum_remaining_read_size) {

            match *preview {
                true => {
                    print_read(read1.name.as_slice(), read1.seq.as_slice(), read1.quals.as_slice(), &base_cuts_read1.1, &base_cuts_read1.2);
                    print_read(read2.name.as_slice(), read2.seq.as_slice(), read2.quals.as_slice(), &base_cuts_read2.1, &base_cuts_read2.2);
                }
                false => {
                    write_read(out_fastq1, &slice_read(&read1, &base_cuts_read1.1, &base_cuts_read1.2)).expect("Unable to write to output file 1.");
                    write_read(out_fastq2, &slice_read(&read2, &base_cuts_read2.1, &base_cuts_read2.2)).expect("Unable to write to output file 2.");
                }
            }
        }
    }
    out_fastq1.flush().expect("Unable to flush output file 1.");
    out_fastq2.flush().expect("Unable to flush output file 2.");
}

fn setup_compressed_file(fastq_output: &Option<String>) -> BufWriter<GzEncoder<Box<dyn Write>>> {
    let writer1: Box<dyn Write> = match fastq_output.clone() {
        Some(file) => Box::new(File::create(file).unwrap()),
        None => Box::new(io::stdout()),
    };

    let out_fastq1 = BufWriter::new(GzEncoder::new(writer1, Compression::default()));
    out_fastq1
}

pub fn write_read(writer: &mut BufWriter<dyn Write>, record: &FastqRecord) -> Result<(), io::Error> {
    writer.write_all(&record.name)?;
    writer.write_all(b"\n")?;
    writer.write_all(&record.seq)?;
    writer.write_all(b"\n+\n")?;
    writer.write_all(&record.quals)?;
    writer.write_all(b"\n")?;
    Ok(())
}

pub fn slice_read(record: &FastqRecord, cut_front: &usize, cut_back: &usize) -> FastqRecord {
    let seq = record.seq.clone();
    let quals = record.quals.clone();
    let name = record.name.clone();
    let seq = &seq[*cut_front..*cut_back];
    let quals = &quals[*cut_front..*cut_back];
    FastqRecord { name, seq: seq.to_vec(), quals: quals.to_vec() }
}


pub fn print_read(name: &[u8], seq: &[u8], qual: &[u8], slice_point_front: &usize, slice_point_rear: &usize) {
    println!(">{}", String::from_utf8(name.to_vec()).unwrap());
    print_format_read(seq, qual, slice_point_front, slice_point_rear);
}

fn print_format_read(seq: &[u8], qual: &[u8], slice_point_front: &usize, slice_point_rear: &usize) {
    let mut seq_string = String::new();
    //seq_string.push_str(&format!("{}",">>>>>>>>>>").white().to_string());
    for index in 0..seq.len() {
        let base = seq.get(index.clone()).unwrap();
        let quality = u32::from(qual[index]);
        assert!(quality >= 33, "Quality score must be at least 33. {} ", quality);

        let mut formated_base = format!("{}", char::from(*base));
        match *base {
            b'A' => { formated_base = formated_base.truecolor(0, color_qual_proportion(&quality), 0).to_string() }
            b'C' => { formated_base = formated_base.truecolor(color_qual_proportion(&quality), 0, 0).to_string() }
            b'G' => { formated_base = formated_base.truecolor(0, 0, color_qual_proportion(&quality)).to_string() }
            b'T' => { formated_base = formated_base.truecolor(color_qual_proportion(&quality), color_qual_proportion(&quality), 0).to_string() }
            _ => {
                let qual_adjusted = color_qual_proportion(&quality);
                formated_base = formated_base.truecolor(qual_adjusted.clone(), qual_adjusted.clone(), qual_adjusted.clone()).to_string()
            }
        }
        if index < *slice_point_front || index >= *slice_point_rear {
            formated_base = formated_base.on_truecolor(255,0,0).to_string();
        } else {
            formated_base = formated_base.on_truecolor(color_qual_proportion(&quality),color_qual_proportion(&quality),color_qual_proportion(&quality)).to_string();
        }

        seq_string.push_str(&format!("{}", formated_base).to_string());
    }

    println!("{}", seq_string);
}

pub fn color_qual_proportion(quality: &u32) -> u8 {
    assert!(quality >= &33, "Quality score must be at least 33. {} ", quality);
    ((((if quality > &93 { 93 } else { quality.clone() }) - 33) as f64 / (93 - 33) as f64) * 255.0) as u8
}