use std::path::PathBuf;
use std::io::{BufRead};

use fastq2comp::io_utils;
use fastq2comp::BaseComp;

#[cfg(tests)]
mod sample_fastq_tests {
    use super::*;
    use test_utils;

    #[test]
    fn test_sample() {
        let reader = test_utils::return_reader(b"@\nAAA\n+\n~~~");
        let mut writer = test_utils::return_writer();

        run(reader, &mut writer, Cli {
            target_read_count: 1u64,
            input: None,
            output: None,
            stdin: false,
            stdout: false,
            header: false,
            seq: false,
            mid: false,
            quals: false,
            min_phred_score: 0,
            n_content: None,
            trimmed_length: Some(2),
        });

        let stringified = writer.get_ref()[0..].to_vec();
        assert_eq!(std::str::from_utf8(&stringified).unwrap().to_string(),
        std::str::from_utf8(b"@\nAA\n+\n~~").unwrap());
    }

    #[test]
    fn test_check_colorspace() {
        let mut read = FASTQRead::new(6);
        let mut reader = return_reader(b"@\nAT1CGN\n+\n!!!!!!");
        read.read(&mut reader);

        assert!(read.check_colorspace())
    }

    #[test]
    fn test_count_n() {
        let mut read = FASTQRead::new(6);
        let mut reader = return_reader(b"@\nNNANNA\n+\n!!!!!!");
        read.read(&mut reader);

        assert_eq!(read.count_n(), 4)
    }

    #[test]
    fn test_get_average_quality() {
        let mut read = FASTQRead::new(6);
        let mut reader = return_reader(b"@\nNNANNA\n+\n!\"#{|}");
        read.read(&mut reader);

        assert_eq!(read.get_average_quality(), 46)
    }
}

use structopt::StructOpt;
#[derive(Debug, StructOpt)]
#[structopt(name = "extract FASTQ base composition", about = "Extracts base composition of FASTQ file and returns result in JSON.")]
struct Cli {
    /// Input file, stdin if not present
    #[structopt(parse(from_os_str), required_unless("stdin"))]
    pub input: Option<PathBuf>,

    /// Output file, use flag -Z to get output to stdout instead
    #[structopt(parse(from_os_str), required_unless("stdout"))]
    pub output: Option<PathBuf>,

    /// Toggles stdin input.
    #[structopt(short = "I", long = "stdin")]
    stdin: bool,

    /// Toggles output to stdout.
    #[structopt(short = "Z", long = "stdout")]
    stdout: bool,

    #[structopt(flatten)]
    sample_args: SampleArgs,
}

#[derive(Debug, StructOpt)]
struct SampleArgs {
    /// Target sample count
    #[structopt()]
    target_read_count: u64,

    /// Sets minimum average quality allowed in sampled reads.
    #[structopt(short = "m", long = "min", default_value="0")]
    min_phred_score: usize,
    /// Sets maximum amount of N's allowed in sample reads.
    #[structopt(short = "n", long = "n-content")]
    n_content: Option<usize>,
    /// Trims each sampled read to given length.
    #[structopt(short = "t", long = "trim")]
    trimmed_length: Option<usize>,
}

use regex::Regex;
/// Abstraction for a single read of FASTQ data
#[derive(Debug)]
struct FASTQRead {
    pub seq: String,
    quals: String,
    // Regex for checking if seq has colorspace
    seq_colorspace_checker: Regex,
}

use log::{debug, info};
impl FASTQRead {
    /// Reads a complete FASTQ statement (composed of 4 lines) into itself
    /// - `reader`: Object implementing `std::io::BufRead` from which to read lines
    /// Note: Will terminate program if EOF reached
    pub fn read<R> (&mut self, reader: &mut R)
    where
        R: BufRead
    {

        //Skips the 1st and 3rd line resp. in 4 lines of input
        for s in [&mut self.seq, &mut self.quals].iter_mut() {
            **s = match reader.lines().skip(1).next() {
                Some(n) => n.expect("Error reading line. Make sure UTF-8 input is supported"),
                None => std::process::exit(0),
            }
        }
    }

    pub fn new (len: usize) -> FASTQRead {
        FASTQRead {
            seq: String::with_capacity(len),
            quals: String::with_capacity(len),

            seq_colorspace_checker: Regex::new(r"\d").unwrap()
        }
    }

    pub fn len (&self) -> usize {
        self.seq.len()
    }

    fn trim(&mut self, len: usize) {
        for s in [&mut self.seq, &mut self.quals].iter_mut() {
            s.truncate(len);
        }
    }

    fn count_n(&self) -> usize {
        self.seq.matches("N").count()
    }

    // Returns true if number is found in seq
    fn check_colorspace(&self) -> bool {
        self.seq_colorspace_checker.is_match(&self.seq)
    }

    fn get_average_quality(&self) -> usize {
        let mut qual_sum: usize = 0;
        for char in self.quals.as_bytes() {
            qual_sum += (*char as usize) - 33;
        }

        qual_sum / self.quals.len()
    }

    fn check_read(read: &mut FASTQRead, args: &SampleArgs) -> bool {
        if let Some(n) = args.trimmed_length {
            read.trim(n);
        };

        // Check for numbers in reads
        if read.check_colorspace() {
            info!("Found numbers in reads - this is probably colorspace");
            std::process::exit(0);
        }

        // Count the N's
        if let Some(n) = args.n_content {
            if read.count_n() > n {
                debug!("N count of current line ({}) too high - skipping", read.count_n());
                return false;
            }
        }

        if read.get_average_quality() < args.min_phred_score {
            debug!("Quality too low - skipping");
            return false;
        }

        true
    }
}

use simple_logger::SimpleLogger;

fn run<R> (args: SampleArgs, mut reader: R) -> String
where R: BufRead {
    let mut read = FASTQRead::new(0);
    read.read(&mut reader);

    // Figure out allotment size based on line size, or provided trim len
    let mut base_comp = BaseComp::init(
        match args.trimmed_length {
            Some(n) => n,
            None => read.len()
        }
    );

    let mut valid_seqs: u64 = 0;

    while valid_seqs <= args.target_read_count {
        if FASTQRead::check_read(&mut read, &args) {
            base_comp.extract(&read.seq);
            valid_seqs += 1;
        }
        read.read(&mut reader);
    }
    info!("Found enough sequences\n");

    for r in base_comp.lib.iter_mut() {
        r.bases.percentage();
    }

    base_comp.jsonify()
}

fn main() {
    //Set up logger.
    SimpleLogger::new().init().unwrap();

    let args = Cli::from_args();
    //Program starts.

    let mut writer = io_utils::get_writer(&args.output);
    let mut reader = io_utils::get_reader(&args.input);

    info!("Arguments recieved: {:#?}", args);

    writeln!(writer, "{}", run(args.sample_args, &mut reader)).expect("Problem printing result");
}