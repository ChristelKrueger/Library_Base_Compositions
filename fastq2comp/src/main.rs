// Reader is a wrapper over BufRead
// Takes in a PathBuf and open it or if no PathBuf is provided, opens up stdin
// And provides an interface over the actual reading.

use std::path::PathBuf;
use std::fs::OpenOptions;
use std::io::{self, BufReader, BufRead, BufWriter, Write};

fn get_reader(input: &Option<PathBuf>) -> Box<dyn BufRead> {
    match input {
        Some(file) => Box::new(BufReader::new(OpenOptions::new().read(true).open(&file).expect("Input file does not exist"))),
        None => {
            let stdin = io::stdin();
            Box::new(BufReader::new(stdin))
        },
    } 
}

use std::io::ErrorKind;
/// Will return writer to File if PathBuf can be opened, will panic if File unavailable
/// And return writer to stdout if PathBuf not given
fn get_writer(output: &Option<PathBuf>) -> BufWriter<Box<dyn Write>> {
    match output {
        Some(file) => BufWriter::new(Box::new(OpenOptions::new().append(true).open(&file).unwrap_or_else(|error| {
            if error.kind() == ErrorKind::NotFound {
                OpenOptions::new().create(true).write(true).open(&file).expect("Problem creating the file!")
            } else {
                panic!("Problem opening the file: {:?}", error);
            }
        }))),
        None => {
            BufWriter::new(Box::new(io::stdout()))
        },
    }
}

use fastq2comp::{run, Cli};
use simple_logger::SimpleLogger;
use structopt::StructOpt;

fn main() {
    //Set up logger.
    SimpleLogger::new().init().unwrap();

    let args = Cli::from_args();
    //Program starts.
    writeln!(get_writer(&args.output), "{}", run(get_reader(&args.input), args)).expect("Problem printing result");
}