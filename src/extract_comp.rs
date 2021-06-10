use std::path::PathBuf;
use std::io::BufRead;

use crate::BaseComp;

#[cfg(test)]
mod sample_fastq_tests {
    use super::*;
    use crate::test_utils::*;

    #[test]
    fn test_json_run() {
        let reader = return_reader(b"@\nAAA\n+\n~~~");
        let args = SampleArgs {
            target_read_count: 1u64,
            min_phred_score: 0,
            n_content: None,
            trimmed_length: Some(2)
        };

        let (result, seqs) = run_json( FASTQReader::new(args, reader));

        assert_eq!(
            result,
            std::str::from_utf8(b"{\"lib\":[{\"pos\":1,\"bases\":{\"A\":100,\"T\":0,\"G\":0,\"C\":0,\"N\":0}},{\"pos\":2,\"bases\":{\"A\":100,\"T\":0,\"G\":0,\"C\":0,\"N\":0}}],\"len\":2}").unwrap()
        );
        assert_eq!(seqs, 1);
    }

    #[test]
    fn test_tsv_run() {
        let reader = return_reader(b"@\nAAA\n+\n~~~");
        let args = SampleArgs {
            target_read_count: 1u64,
            min_phred_score: 0,
            n_content: None,
            trimmed_length: Some(2)
        };

        let (result, seqs) = run_tsv( FASTQReader::new(args, reader));

        assert_eq!(
            result,
            std::str::from_utf8(b"100\t0\t0\t0\t0\t100\t0\t0\t0\t0").unwrap()
        );
        assert_eq!(seqs, 1);
    }

    #[test]
    fn test_check_colorspace() {
        let mut read = FASTQRead::new(6);
        let mut reader = return_reader(b"@\nAT1CGN\n+\n!!!!!!");
        read.read(&mut reader, None);

        assert!(read.check_colorspace("AT1CGN"))
    }

    #[test]
    fn test_count_n() {
        let mut read = FASTQRead::new(6);
        let mut reader = return_reader(b"@\n\n+\n!!!!!!");
        read.read(&mut reader, None);

        assert_eq!(FASTQRead::count_n("NNANNA"), 4)
    }

    #[test]
    fn test_get_average_quality() {
        assert_eq!(FASTQRead::get_average_quality("#{|}"), 68)
    }
}

use structopt::StructOpt;
#[derive(Debug, StructOpt)]
#[structopt(name = "extract FASTQ base composition", about = "Extracts base composition of FASTQ file and returns result in JSON.")]
pub struct Cli {
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

#[derive(Debug, StructOpt)]
pub struct SampleArgs {
    /// Target sample count
    #[structopt()]
    pub target_read_count: u64,

    /// Sets minimum average quality allowed in sampled reads.
    #[structopt(short = "m", long = "min", default_value="0")]
    pub min_phred_score: usize,
    /// Sets maximum amount of N's allowed in sample reads.
    #[structopt(short = "n", long = "n-content")]
    pub n_content: Option<usize>,
    /// Trims each sampled read to given length.
    #[structopt(short = "t", long = "trim")]
    pub trimmed_length: Option<usize>,
}

use regex::Regex;

use lazy_static::lazy_static;
lazy_static! {
    // Regex for checking if seq has numbers
    static ref SEQCOLORSPACECHECKER: Regex = Regex::new(r"\d").unwrap();
}

/// Abstraction for a single read of FASTQ data
#[derive(Debug)]
pub(crate) struct FASTQRead {
    pub seq: String,
    quals: String,
}

use log::{debug, info};
impl FASTQRead {

    /// Reads a complete FASTQ statement (composed of 4 lines) into itself
    /// - `reader`: Object implementing `std::io::BufRead` from which to read lines
    /// - Returns `None` if EOF reached.
    fn read_fastq (&mut self, reader: &mut impl BufRead) -> Option<()> {
        //Skips the 1st and 3rd line resp. in 4 lines of input
        for s in [&mut self.seq, &mut self.quals].iter_mut() {
            **s = match reader.lines().nth(1) {
                Some(n) => n.expect("Error reading line. Make sure UTF-8 input is supported"),
                None => {debug!("Input reading finished"); return None},
            }
        }

        Some(())
    }

    /// Reads a complete FASTQ statement (composed of 4 lines) into itself.
    /// - *Important*: Panics if read length of strings is not same as strings read before.
    /// - `reader`: Object implementing `std::io::BufRead` from which to read lines
    /// - Returns `None` if EOF reached.
    fn read (&mut self, reader: &mut impl BufRead, len: Option<usize>) -> Option<()>
    {
        let (seq_len, quals_len) = (self.seq.len(), self.quals.len());
        
        let result = self.read_fastq(reader);

        // If there have been no reads, then perform first read without checking for lengths
        if seq_len == 0 || quals_len == 0 {
            return result;
        }

        // behaviour is meant to panic only when uneven reads cannot be fixed by trimming
        if len.unwrap_or(usize::MAX) > seq_len && (seq_len != self.seq.len() || quals_len != self.quals.len()) {
            panic!("Reads have inconsistent lengths. Particularly this read: \n{}\n{}\n \
            Expected read to be: {} length but it was: {} length", self.seq, self.quals, seq_len, self.seq.len());
        }

        result
    }

    fn new (len: usize) -> FASTQRead {
        FASTQRead {
            seq: String::with_capacity(len),
            quals: String::with_capacity(len),
        }
    }

    fn count_n(seq: &str) -> usize {
        seq.matches('N').count()
    }

    // Returns true if number is found in seq
    fn check_colorspace(&self, seq: &str) -> bool {
        SEQCOLORSPACECHECKER.is_match(seq)
    }

    fn get_average_quality(quals: &str) -> usize {
        let mut qual_sum: usize = 0;
        for char in quals.as_bytes() {
            qual_sum += (*char as usize) - 33;
        }

        qual_sum / quals.len()
    }

    fn trim(str: &str, len: Option<usize>) -> &str {
        match len {
            Some(n) => &str[0..n],
            None => &str[0..],
        }
    }

    fn check_read(&mut self, args: &SampleArgs) -> bool {
        let seq = FASTQRead::trim(&self.seq, args.trimmed_length);
        let quals = FASTQRead::trim(&self.seq, args.trimmed_length);

        // Check for numbers in reads
        if self.check_colorspace(seq) {
            info!("Found numbers in reads - this is probably colorspace");
            std::process::exit(0);
        }

        // Count the N's
        if let Some(n) = args.n_content {
            if FASTQRead::count_n(seq) > n {
                debug!("N count of current line ({}) too high - skipping", FASTQRead::count_n(seq));
                return false;
            }
        }

        if FASTQRead::get_average_quality(quals) < args.min_phred_score {
            debug!("Quality too low - skipping");
            return false;
        }

        true
    }
}

use reservoir_sampling::unweighted::l as sample;

/// Takes in reader (for FASTQ lines) and SampleArgs, returns JSONified string and
/// total number of reads processed after applying SampleArgs.
pub fn run_json<T> (fastq_reader: FASTQReader<T>) -> (String, u64)
where T: BufRead
{
    let (mut comp, lines_read) = run_core (fastq_reader);

    (comp.jsonify(), lines_read)
}

pub fn run_tsv<T> (fastq_reader: FASTQReader<T>) -> (String, u64)
where T: BufRead
{
    let (comp, lines_read) = run_core (fastq_reader);

    ({let mut s = comp.lib.into_iter().flat_map(|b| b.bases.iter()).
        fold(String::new(), |acc, curr| acc + &curr.to_string() + "\t");
        s.pop(); // remove trailing ',' to make it valid tsv
        s
    },
    lines_read)
}

/// Takes in reader (for FASTQ lines) and SampleArgs, returns [`BaseComp`] and
/// total number of reads processed after applying SampleArgs.
fn run_core<T> (fastq_reader: FASTQReader<T>) -> (BaseComp, u64)
where T: BufRead
{
    //TODO: Convert args.target_read_count to usize or figure out how to allocate u64-sized vec
    let sampled_seqs = fastq_reader.sample_random();

    // Figure out allotment size based on line size, or provided trim len
    let mut base_comp = BaseComp::init(sampled_seqs[0].len());

    let mut lines_read = 0u64;
    for seq in sampled_seqs {
        if seq.is_empty() {
            break;
        }
        base_comp.extract(&seq);
        lines_read += 1;
    }

    for r in base_comp.lib.iter_mut() {
        r.bases.percentage();
    }

    (base_comp, lines_read)
}


pub struct FASTQReader<T: BufRead> {
    curr: FASTQRead,
    reader: T,
    sample_args: SampleArgs,
    pub curr_valid_reads: u64,
    pub target_read_count: u64,
}

impl<T: BufRead> FASTQReader<T> {
    pub fn new (args: SampleArgs, reader: T) -> FASTQReader<T> {
        let read = FASTQRead::new(args.trimmed_length.unwrap_or(0));
        let target_read_count = args.target_read_count.clone();

        FASTQReader {
            curr: read,
            reader,
            sample_args: args,
            curr_valid_reads: 0,
            target_read_count,
        }
    }
    pub fn sample_random (self) -> Vec<String> {
        let mut sampled_seqs = vec![String::new(); self.target_read_count as usize];

        // Randomly sample FASTQ reads
        sample(self, sampled_seqs.as_mut_slice());
        sampled_seqs
    }
}

impl<T: BufRead> Iterator for FASTQReader<T> {
    type Item = String;

    fn next (&mut self) -> Option<String> {
        if self.curr_valid_reads >= self.sample_args.target_read_count {
            return None
        }

        self.curr.read(&mut self.reader, self.sample_args.trimmed_length);
        while !FASTQRead::check_read(&mut self.curr, &self.sample_args) {
            if self.curr.read(&mut self.reader, self.sample_args.trimmed_length).is_none() {
                return None;
            }
        }
        self.curr_valid_reads += 1;

        Some(FASTQRead::trim(&self.curr.seq, self.sample_args.trimmed_length).to_string())
    }

}
