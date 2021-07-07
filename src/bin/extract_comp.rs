use fastq2comp::extract_comp::{Cli, FASTQReader, run_tsv, run_json};
use fastq2comp::io_utils;

use structopt::StructOpt;

fn main() {
    let args = Cli::from_args();
    //Program starts.

    let mut writer = io_utils::get_writer(&args.output);
    let mut reader = io_utils::get_reader(&args.input, args.compressed);

    let (result, reads) = 
        if args.tsv {run_tsv(FASTQReader::new(args.sample_args, &mut reader))}
        else {run_json(FASTQReader::new(args.sample_args, &mut reader))};

    writeln!(writer, "{}", result).expect("Problem printing result");
    eprintln!("Extracted base composition of {} reads.", reads);
}