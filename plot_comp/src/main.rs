use plot_comp::{run, Cli};
use fastq_io;
use simple_logger::SimpleLogger;
use structopt::StructOpt;

fn main() {
    //Set up logger.
    SimpleLogger::new().init().unwrap();

    let args = Cli::from_args();

    run(fastq_io::get_reader(&args.input), args.libs).expect("Error drawing chart");
}