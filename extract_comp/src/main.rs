use extract_comp::{run, Cli};
use fastq_io;
use simple_logger::SimpleLogger;
use structopt::StructOpt;
use std::io::{self, BufWriter};

fn main() {
    //Set up logger.
    SimpleLogger::new().init().unwrap();

    let args = Cli::from_args();
    //Program starts.
    run(fastq_io::get_reader(&args.input), BufWriter::new(io::stdout()));
}