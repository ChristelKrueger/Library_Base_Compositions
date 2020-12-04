use fastq2comp::extract_fastq_comp::{SampleArgs, run};
use fastq2comp::io_utils;
use std::path::PathBuf;

/// This example extracts the base composition of a file
/// and prints it JSON format.

fn main() {
    let mut reader = io_utils::get_reader(&Some(PathBuf::from("sample_fastq_file.fastq")));
    println!{"{}",
        run(SampleArgs {
            trimmed_length: Some(50),
            target_read_count: 10,
            min_phred_score: 0,
            n_content: None
        }, &mut reader)
    };
}