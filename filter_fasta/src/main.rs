use std::env;
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, Write};
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 5 {
        eprintln!(
            "Usage: {} -i <input.fasta> -l <list.txt> [-o <output.fasta>]",
            args[0]
        );
        process::exit(1);
    }

    let input_fasta = &args[2];
    let list_file = &args[4];
    let output_fasta = if args.len() > 5 {
        &args[6]
    } else {
        "filtered.fasta"
    };

    // Open the input and output files
    let in_fasta_file = match File::open(input_fasta) {
        Ok(file) => file,
        Err(_) => {
            eprintln!("Error opening input FASTA file: {}", input_fasta);
            process::exit(1);
        }
    };
    let out_fasta_file = match OpenOptions::new()
        .write(true)
        .create(true)
        .open(output_fasta)
    {
        Ok(file) => file,
        Err(_) => {
            eprintln!("Error opening/creating output FASTA file: {}", output_fasta);
            process::exit(1);
        }
    };

    let sample_list = match File::open(list_file) {
        Ok(file) => file,
        Err(_) => {
            eprintln!("Error opening sample list file: {}", list_file);
            process::exit(1);
        }
    };
    let sample_list_reader = BufReader::new(sample_list);
    let samples_to_keep: Vec<String> = sample_list_reader.lines().filter_map(Result::ok).collect();

    // Read the input FASTA
    let in_fasta_reader = BufReader::new(in_fasta_file);
    let in_fasta_lines: Vec<String> = in_fasta_reader.lines().filter_map(Result::ok).collect();

    // Create a buffered writer for the output file
    let mut out_fasta_writer = std::io::BufWriter::new(out_fasta_file);

    // Process the FASTA file and filter records sequentially
    let mut keep_record = false;
    let mut seq_buffer = String::new();
    let mut header: Option<String> = None;

    // Iterate over lines and handle multi-line FASTA records
    for line in in_fasta_lines {
        if line.starts_with('>') {
            // Write the previous record if it was selected
            if let Some(ref h) = header {
                if keep_record {
                    let _ = writeln!(out_fasta_writer, "{}", h); // Write header
                    let _ = writeln!(out_fasta_writer, "{}", seq_buffer); // Write sequence
                }
            }
            // Start a new record
            header = Some(line.clone());
            seq_buffer.clear();
            let record_name: &String = &String::from(&line[1..]);
            keep_record = samples_to_keep.contains(record_name);
        } else if keep_record {
            // Accumulate sequence lines for multi-line records
            seq_buffer.push_str(&line);
        }
    }

    // Handle the last record if it's kept
    if let Some(ref h) = header {
        if keep_record {
            let _ = writeln!(out_fasta_writer, "{}", h); // Write header
            let _ = writeln!(out_fasta_writer, "{}", seq_buffer); // Write sequence
        }
    }

    println!("Filtered FASTA written to: {}", output_fasta);
}
