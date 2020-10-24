use sample_fastq_file::sample::{run, Cli};
use sample_fastq_file::fastq_io;
use simple_logger::SimpleLogger;
use structopt::StructOpt;

fn main() {
    //Set up logger.
    SimpleLogger::new().init().unwrap();

    let args = Cli::from_args();
    //Program starts.
    run(fastq_io::get_reader(&args.input), fastq_io::get_writer(&args.output), args);
}