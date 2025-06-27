use std::{
    fs::{self, File},
    io::BufReader,
    path::{Path, PathBuf},
    str,
};

use clap::Parser as ClapParser;
use fastq::{each_zipped, Parser, Record};
use flate2::read::MultiGzDecoder;
use indicatif::{ProgressBar, ProgressStyle};
use polars::prelude::ParquetWriter;
use polars::prelude::*;
use rayon::prelude::*;

#[derive(Debug)]
struct ReadInfo {
    read_id: String,
    primer_fwd: String,
    primer_rev: String,
    r1_length: usize,
    r2_length: usize,
    umi5: String,
    umi3: String,
}

/// CLI for processing FASTQ files into Parquet summaries
#[derive(ClapParser, Debug)]
#[command(author, version, about)]
struct Args {
    /// Input directory containing *_1.fq.gz and *_2.fq.gz
    #[arg(short, long)]
    input: String,

    /// Output directory to write .parquet files
    #[arg(short, long)]
    output: String,
}

fn process_sample(read1_path: &Path, output_dir: &Path) {
    let sample = read1_path
        .file_name()
        .unwrap()
        .to_string_lossy()
        .split('_')
        .next()
        .unwrap()
        .to_string();

    let output_file = output_dir.join(format!("{sample}.parquet"));
    if output_file.exists() {
        println!("Skipping {}: already processed", output_file.display());
        return;
    }

    let read2_path = PathBuf::from(read1_path.to_string_lossy().replace("_1.fq.gz", "_2.fq.gz"));
    let r1_reader = MultiGzDecoder::new(File::open(read1_path).unwrap());
    let r2_reader = MultiGzDecoder::new(File::open(&read2_path).unwrap());

    let parser1 = Parser::new(BufReader::new(r1_reader));
    let parser2 = Parser::new(BufReader::new(r2_reader));

    let mut records = Vec::new();

    each_zipped(parser1, parser2, |rec1_opt, rec2_opt| {
        match (rec1_opt, rec2_opt) {
            (Some(rec1), Some(rec2)) => {
                // Print human-readable header from rec1
                let head1_str = str::from_utf8(rec1.head()).unwrap_or("[invalid utf8]");
                // Parse ID and desc
                let mut head_parts = head1_str.split_whitespace();
                let id1 = head_parts.next().unwrap_or("[missing_id]");
                let desc = head_parts.next().unwrap_or("[missing_desc]");

                let head2_str = str::from_utf8(rec2.head()).unwrap_or("[invalid utf8]");
                let mut head2_parts = head2_str.split_whitespace();
                let id2 = head2_parts.next().unwrap_or("[missing_id]");

                if id1 != id2 {
                    eprintln!("Mismatched IDs: {} vs {}", id1, id2);
                    return (false, false);
                }

                // Parse primer info
                let (primer_fwd, primer_rev) = match desc.split_once('+') {
                    Some(pair) => pair,
                    None => {
                        eprintln!("Invalid desc format: {}", desc);
                        return (false, false);
                    }
                };

                // Parse UMI5/UMI3 from read ID
                let mut id_parts = id1.split('_');
                let _read_id = id_parts.next().unwrap_or("");
                let umi5 = id_parts.next().unwrap_or("").to_string();
                let umi3 = id_parts.next().unwrap_or("").to_string();

                records.push(ReadInfo {
                    read_id: id1.to_string(),
                    primer_fwd: primer_fwd.to_string(),
                    primer_rev: primer_rev.to_string(),
                    r1_length: rec1.seq().len(),
                    r2_length: rec2.seq().len(),
                    umi5,
                    umi3,
                });

                (true, true)
            }

            (None, None) => {
                // EOF on both files: clean finish
                (false, false)
            }

            _ => {
                // One file ended early
                eprintln!("Error: FASTQ files are not aligned — one file has more reads");
                (false, false)
            }
        }
    })
    .expect("FASTQ parsing failed");

    if !records.is_empty() {
        let df = df![
            "read_id" => records.iter().map(|r| r.read_id.as_str()).collect::<Vec<_>>(),
            "primer_fwd" => records.iter().map(|r| r.primer_fwd.as_str()).collect::<Vec<_>>(),
            "primer_rev" => records.iter().map(|r| r.primer_rev.as_str()).collect::<Vec<_>>(),
            "r1_length" => records.iter().map(|r| r.r1_length as u32).collect::<Vec<_>>(),
            "r2_length" => records.iter().map(|r| r.r2_length as u32).collect::<Vec<_>>(),
            "umi5" => records.iter().map(|r| r.umi5.as_str()).collect::<Vec<_>>(),
            "umi3" => records.iter().map(|r| r.umi3.as_str()).collect::<Vec<_>>(),
        ]
        .unwrap();

        let mut df = df;
        let file = File::create(output_file).unwrap();
        ParquetWriter::new(file).finish(&mut df).unwrap();
    }
}

fn main() {
    let args = Args::parse();

    let data_dir = PathBuf::from(&args.input);
    let output_dir = PathBuf::from(&args.output);
    fs::create_dir_all(&output_dir).unwrap();

    let mut read1_files: Vec<_> = fs::read_dir(&data_dir)
        .unwrap()
        .filter_map(|e| {
            let path = e.unwrap().path();
            if path.extension()? == "gz" && path.file_name()?.to_string_lossy().contains("_1.fq.gz")
            {
                Some(path)
            } else {
                None
            }
        })
        .collect();
    read1_files.sort_unstable_by_key(|p| p.file_name().map(|s| s.to_os_string()));

    rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build_global()
        .unwrap();

    let pb = ProgressBar::new(read1_files.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
            .unwrap(),
    );

    read1_files.par_iter().for_each(|read1_path| {
        process_sample(read1_path, &output_dir);
        pb.inc(1);
    });

    pb.finish_with_message("Done");
}
