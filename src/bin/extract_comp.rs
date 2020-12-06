use fastq2comp::extract_comp::{Cli, run};
use fastq2comp::io_utils;

use log::info;
use simple_logger::SimpleLogger;
use structopt::StructOpt;

fn main() {
    //Set up logger.
    SimpleLogger::new().init().unwrap();

    let args = Cli::from_args();
    //Program starts.

    let mut writer = io_utils::get_writer(&args.output);
    let mut reader = io_utils::get_reader(&args.input);

    info!("Arguments recieved: {:#?}", args);

    let (result, seqs) = run(args.sample_args, &mut reader);

    writeln!(writer, "{}", result).expect("Problem printing result");
    eprintln!("Extracted base composition of {} reads.", seqs);
}