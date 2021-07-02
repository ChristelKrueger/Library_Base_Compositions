use fastq2comp::extract_comp::{SampleArgs, run_json};
use fastq2comp::io_utils;
use fastq2comp::extract_comp::FASTQReader;

use std::{path::PathBuf, io::Write};
use std::fs::File;

/// This example extracts the base composition of a file
/// and prints it JSON format.

fn main() {
    let mut reader = io_utils::get_reader(&Some(PathBuf::from("examples/extract-comp/in.fastq")), false);

    let mut file = match File::create(&PathBuf::from("examples/extract-comp/out.json")) {
        Err(why) => panic!("Couldn't open output JSON file: {}", why),
        Ok(file) => file,
    };

    let result = run_json(FASTQReader::new( SampleArgs {
        trimmed_length: Some(50),
        target_read_count: 10,
        min_phred_score: 0,
        n_content: None,
    }, &mut reader));

    match file.write_all(result.0.as_bytes()) {
        Err(why) => panic!("couldn't write to output JSON file: {}", why),
        Ok(_) => println!("successfully wrote to output JSON file, read {} reads", result.1),
    }
}