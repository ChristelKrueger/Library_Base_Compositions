use fastq2comp::extract_comp::{SampleArgs, run};
use fastq2comp::io_utils;
use std::{path::PathBuf, io::Write};
use std::fs::File;

/// This example extracts the base composition of a file
/// and prints it JSON format.

fn main() {
    let mut reader = io_utils::get_reader(&Some(PathBuf::from("examples/extract-comp/in.fastq")));

    let mut file = match File::create(&PathBuf::from("examples/extract-comp/out.json")) {
        Err(why) => panic!("Couldn't open output JSON file: {}", why),
        Ok(file) => file,
    };

    let result = run(SampleArgs {
        trimmed_length: Some(50),
        target_read_count: 10,
        min_phred_score: 0,
        n_content: None
    }, &mut reader).0;

    match file.write_all(result.as_bytes()) {
        Err(why) => panic!("couldn't write to output JSON file: {}", why),
        Ok(_) => println!("successfully wrote to output JSON file"),
    }
}