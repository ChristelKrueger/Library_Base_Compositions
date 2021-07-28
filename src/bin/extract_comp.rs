use std::path::PathBuf;

use fastq2comp::extract_comp::{FASTQReader, run_tsv, run_json, SampleArgs};
use fastq2comp::io_utils;

use structopt::StructOpt;

#[derive(Debug, StructOpt)]
#[structopt(name = "extract FASTQ base composition", about = "Extracts base composition of FASTQ file and returns result in JSON.")]
struct Cli {
    /// Toggles stdin input.
    #[structopt(short = "I", long = "stdin")]
    stdin: bool,

    /// Toggles output to stdout.
    #[structopt(short = "Z", long = "stdout")]
    stdout: bool,

    /// Toggles output in tsv
    #[structopt(long = "tsv")]
    pub tsv: bool,

    #[structopt(flatten)]
    pub sample_args: SampleArgs,

    /// Input file, use flag -I to get output to stdout instead
    #[structopt(parse(from_os_str), required_unless("stdin"))]
    pub input: Option<PathBuf>,

    /// Output file, use flag -Z to get output to stdout instead
    #[structopt(parse(from_os_str), required_unless("stdout"))]
    pub output: Option<PathBuf>,

    /// Uncompress file using gzip
    #[structopt(short = "C", long = "gzip")]
    pub compressed: bool,
}

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