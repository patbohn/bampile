extern crate bam;
extern crate bio;
extern crate clap;

use bam::{RecordReader, RecordWriter};
use bio::io::fasta;
use clap::{Arg, ArgAction, Command};
use std::collections::{HashMap, BTreeSet};
use std::fs::{File, create_dir_all};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let matches = Command::new("BAM Pileup Analyzer")
        .version("1.0")
        .author("Patrick Bohn")
        .about("Counts matches/mismatches at specified positions in a BAM file")
        .arg(
            Arg::new("bam")
                .short('b')
                .long("bam")
                .value_name("BAM_FILE")
                .help("Input BAM file")
                .required(true),
        )
        .arg(
            Arg::new("bed")
                .short('e')
                .long("bed")
                .value_name("BED_FILE")
                .help("Input BED file")
                .required(true),
        )
        .arg(
            Arg::new("fasta")
                .short('f')
                .long("fasta")
                .value_name("FASTA_FILE")
                .help("Reference FASTA file")
                .required(true),
        )
        .arg(
            Arg::new("output_dir")
                .short('o')
                .long("output-dir")
                .value_name("OUTPUT_DIR")
                .help("Output directory for TSV.gz files")
                .required(true),
        )
        .arg(
            Arg::new("qscore_cutoff")
                .short('q')
                .long("qscore")
                .value_name("QSCORE")
                .help("Minimum Q-score cutoff for a match")
                .default_value("30"),
        )
        .get_matches();

    let bam_file_path = matches.get_one::<String>("bam").unwrap();
    let bed_file_path = matches.get_one::<String>("bed").unwrap();
    let fasta_file_path = matches.get_one::<String>("fasta").unwrap();
    let output_dir_path = matches.get_one::<String>("output_dir").unwrap();
    let qscore_cutoff: u8 = matches
        .get_one::<String>("qscore_cutoff").unwrap()
        .parse()
        .expect("Invalid Q-score cutoff");

    // Create the output directory if it doesn't exist
    create_dir_all(output_dir_path)?;

    // Open the BAM file for reading
    let mut bam = bam::IndexedReader::from_path(&Path::new(bam_file_path)).unwrap();

    // Open the reference genome FASTA file
    let reference = fasta::Reader::from_path(fasta_file_path)?;

    // Load the regions of interest from the BED file into a HashMap
    let regions_of_interest = load_bed_regions(bed_file_path)?;

    // Create a HashMap to store counts for each reference sequence and read_id
    let mut reference_counts: HashMap<String, HashMap<String, (usize, usize)>> = HashMap::new();

    // Iterate through regions of interest and extract pileup
    for (ref_name, (start, end)) in regions_of_interest {
        let ref_id = bam.header().reference_id(&ref_name).unwrap();

        let pileup = bam.fetch(&bam::Region::new(ref_id, start, end)).unwrap();

        let mut reference_sequence = Vec::new();
        reference.fetch_into_vec(ref_name, start as u64, (end - 1) as u64, &mut reference_sequence)?;
        
        
        
        let mut record = bam::Record::new();
        loop {
           // reader: impl RecordReader
           // New record is saved into record.
            match pileup.read_into(&mut record) {
               // No more records to read.
                Ok(false) => break,
                Ok(true) => {},
                Err(e) => panic!("{}", e),
            }
            // Do somethind with the record.
            let record = record?;
            let seq = record.seq();
            let qual = record.qual();
            let read_id = record.read().qname();
            let reference_base = reference_sequence[(record.pos() - start) as usize];

            let (num_matches, num_mismatches) = count_matches_mismatches(seq, qual, reference_base, qscore_cutoff);
            let ref_counts = reference_counts
            .entry(ref_name.to_string())
            .or_insert_with(HashMap::new);

            let (matches, mismatches) = ref_counts
                .entry(read_id.to_string())
                .or_insert((0, 0));
            *matches += num_matches;
            *mismatches += num_mismatches;
        }
    }

    // Write results to summary files, one per reference sequence
    for (ref_name, read_counts) in reference_counts {
        let output_file_name = format!("{}/{}.tsv.gz", output_dir_path, sanitize_filename(&ref_name));
        let output_file = File::create(output_file_name)?;
        let mut output_writer =
            flate2::write::GzEncoder::new(output_file, flate2::Compression::default());
        writeln!(output_writer, "read_id\tnum_matches\tnum_mismatches")?;
        for (read_id, (matches, mismatches)) in read_counts {
            writeln!(output_writer, "{}\t{}\t{}", read_id, matches, mismatches)?;
        }
    }

    Ok(())
}

fn count_matches_mismatches(
    sequence: &[u8],
    qual: &[u8],
    reference_base: u8,
    qscore_cutoff: u8,
) -> (usize, usize) {
    let mut num_matches = 0;
    let mut num_mismatches = 0;

    for (base, qscore) in sequence.iter().zip(qual.iter()) {
        if *base == reference_base || *qscore >= qscore_cutoff {
            num_matches += 1;
        } else {
            num_mismatches += 1;
        }
    }

    (num_matches, num_mismatches)
}

fn load_bed_regions(bed_file_path: &str) -> Result<HashMap<String, (u32, u32)>, Box<dyn std::error::Error>> {
    let bed_file = File::open(bed_file_path)?;
    let reader = BufReader::new(bed_file);
    let mut regions_of_interest: HashMap<String, (u32, u32)> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.trim().split('\t').collect();

        if fields.len() >= 3 {
            let chromosome = fields[0].to_string();
            let start = fields[1].parse::<u32>()?;
            let end = fields[2].parse::<u32>()?;
            regions_of_interest.insert(chromosome, (start, end));
        }
    }

    Ok(regions_of_interest)
}

fn sanitize_filename(filename: &str) -> String {
    filename.chars().filter(|c| c.is_alphanumeric() || *c == '_' || *c == '-').collect()
}
