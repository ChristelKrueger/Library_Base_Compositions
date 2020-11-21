#[cfg(tests)]
mod tests {
    use super::*;
    use fastq_io::{return_reader, return_writer};

    #[test]
    fn test_sample() {
        let reader = fastq_io::test_utils::return_reader(b"@\nAAA\n+\n~~~");
        let mut writer = fastq_io::test_utils::return_writer();

        sample::run(reader, &mut writer, sample::Cli {
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

use log::{debug, info};
use structopt::StructOpt;
use std::path::PathBuf;
use std::io::{BufRead, Write};

#[derive(Debug, StructOpt)]
#[structopt(name = "sample fastq file", about = "Filters given FASTQ file according to given specifications.")]
pub struct Cli {
    /// Target sample count
    #[structopt()]
    target_read_count: u64,

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
    /// Don't include header in output
    #[structopt(long = "skip-header")]
    header: bool,

    /// Don't include seq in output
    #[structopt(long = "skip-seq")]
    seq: bool,

    /// Don't include mid in output
    #[structopt(long = "skip-mid")]
    mid: bool,

    /// Don't include quals in output
    #[structopt(long = "skip-quals")]
    quals: bool,

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
// Abstraction for a single read of FASTQ data
#[derive(Debug)]
pub struct FASTQRead {
    pub header: String,
    pub seq: String,
    pub mid: String,
    pub quals: String,
    // Regex for checking if seq has colorspace
    seq_colorspace_checker: Regex,
}

impl FASTQRead {
    /// Reads a complete FASTQ statement (composed of 4 lines) into itself
    /// - `reader`: Object implementing `std::io::BufRead` from which to read lines
    /// Note: Will terminate program if EOF reached
    pub fn read<R> (&mut self, reader: &mut R)
    where
        R: BufRead
    {
        for s in [&mut self.header, &mut self.seq, &mut self.mid, &mut self.quals].iter_mut() {
            s.clear();
            if reader.read_line(*s).expect("Error reading line. Make sure terminal supports UTF-8 input") == 0 {
                process::exit(0);
            }

            //Trims string to eliminate newline.
            **s = s.trim().to_string();
        }
    }
    
    /// Prints the contents of FASTQRead. Pass false for a flag if it shouldn't be printed
    /// - `writer` Object implementing `std::io::Write` to which to write FASTQ-formatted data
    /// - `header` Pass true to print the _header_ line of a FASTQ Read (starts with `@`)
    /// - `seq` Pass true to print the _sequence_ line of a FASTQ Read
    /// - `mid` Pass true to print the _mid_ line of a FASTQ Read (starts with `+`)
    /// - `quals` Pass true to print the _quality_ line of a FASTQ Read
    pub fn write<W> (&self, writer: &mut W, header: bool, seq: bool, mid: bool, quals: bool)
    where
        W: Write,
    {
        if header {write!(writer, "{}\n", self.header).expect("Error writing to file/stdout.");}
        if seq {write!(writer, "{}\n", self.seq).expect("Error writing to file/stdout.");}
        if mid {write!(writer, "{}\n", self.mid).expect("Error writing to file/stdout.");}
        if quals {write!(writer, "{}\n", self.quals).expect("Error writing to file/stdout.");}
    }

    pub fn new (len: usize) -> FASTQRead {
        FASTQRead {
            header:  String::with_capacity(len),
            seq: String::with_capacity(len),
            mid: String::with_capacity(len),
            quals: String::with_capacity(len),

            seq_colorspace_checker: Regex::new(r"\d").unwrap()
        }
    }

    pub fn trim(&mut self, len: usize) {
        for s in [&mut self.seq, &mut self.quals].iter_mut() {
            s.truncate(len);
        }
    }

    pub fn count_n(&self) -> usize {
        self.seq.matches("N").count()
    }

    // Returns true if number is found in seq
    pub fn check_colorspace(&self) -> bool {
        self.seq_colorspace_checker.is_match(&self.seq)
    }

    pub fn get_average_quality(&self) -> usize {
        let mut qual_sum: usize = 0;
        for char in self.quals.as_bytes() {
            qual_sum += (*char as usize) - 33;
        }

        qual_sum / self.quals.len()
    }
}

use std::process;
fn check_read(read: &mut FASTQRead, args: &SampleArgs) -> bool {
    if let Some(n) = args.trimmed_length {
        read.trim(n);
    };

    // Check for numbers in reads
    if read.check_colorspace() {
        info!("Found numbers in reads - this is probably colorspace");
        process::exit(0);
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

pub fn run<R, W>(mut reader: R, mut writer: W, args: Cli)
where
    R: BufRead,
    W: Write,
{
    info!("Arguments recieved: {:#?}", args);

    let mut read = FASTQRead::new(0);

    let mut valid_seqs: u64 = 0;

    // Set up write buffer to speed up write performance.
    // As we won't have to flush everytime we write.
    // TODO: There is an issue with BufWriter, it is not flushing automatically,
    // and flushing at the end does not flush the data in buffer inserted during the loop.
    // So as a stopgap we are flushing after printing a single read, which defeats the purpose of
    // using BufWriter

    while valid_seqs <= args.target_read_count {
        read.read(&mut reader);
        if check_read(&mut read, &args.sample_args) {
            read.write(&mut writer, !args.sample_args.header, !args.sample_args.seq, !args.sample_args.mid, !args.sample_args.quals);
            writer.flush().expect("Error flushing stream");

            valid_seqs += 1;
        }
    }
    info!("Found enough sequences\n");
}