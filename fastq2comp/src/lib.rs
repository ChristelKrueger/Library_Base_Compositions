#[cfg(tests)]
mod sample_fastq_tests {
    use super::*;

    use std::io::Cursor;
    /// Returns reader which implements Read trait.
    /// - `s`: Data which should be yielded by the reader upon read
    fn return_reader<'a> (s: &'a[u8]) -> Cursor<&[u8]> {
        Cursor::new(s)
    }

    /// Returns a writer which implements Write trait.
    /// `writer.get_ref()[0..].to_vec()` can be used to get the data written to the writer.
    fn return_writer() -> Cursor<Vec<u8>> {
        Cursor::new(Vec::<u8>::new())
    }

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
use std::io::BufRead;


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
    pub seq: String,
    quals: String,
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

pub fn run<R>(mut reader: R, args: Cli) -> String
where
    R: BufRead,
{
    info!("Arguments recieved: {:#?}", args);

    let mut read = FASTQRead::new(0);
    read.read(&mut reader);

    // Figure out allotment size based on line size, or provided trim len
    let mut base_comp = BaseComp::init(
        match args.sample_args.trimmed_length {
            Some(n) => n,
            None => read.len()
        }
    );

    let mut valid_seqs: u64 = 0;

    while valid_seqs <= args.target_read_count {
        if FASTQRead::check_read(&mut read, &args.sample_args) {
            base_comp.extract(&read.seq);
            valid_seqs += 1;
        }
        read.read(&mut reader);
    }
    info!("Found enough sequences\n");

    for r in base_comp.lib.iter_mut() {
        r.percentage();
    }

    base_comp.jsonify()
}

use serde::{Serialize, Deserialize};
use serde_json;

#[derive(Serialize, Deserialize)]
#[allow(non_snake_case)]
pub struct ColBaseComp {
    pub pos: usize,
    pub A: usize,
    pub T: usize,
    pub G: usize,
    pub C: usize,
    pub N: usize,
}

#[cfg(test)]
mod col_base_comp_tests {
    use super::*;
    #[test]
    fn test_extract() {
        let mut read = ColBaseComp::new(0);
        read.extract(&b'A');
        assert_eq!(read.A, 1);

        read.extract(&b'C');
        assert_eq!(read.C, 1);

        read.extract(&b'T');
        assert_eq!(read.T, 1);

        read.extract(&b'G');
        assert_eq!(read.G, 1);

        read.extract(&b'N');
        assert_eq!(read.N, 1);
    }
    #[test]
    fn test_percentage() {
        let mut read = ColBaseComp::new(0);
        for c in "AGTCGA".as_bytes().iter() {
            read.extract(c);
        }
        read.percentage();

        assert_eq!(read.A, 20);
        assert_eq!(read.C, 20);
        assert_eq!(read.T, 20);
        assert_eq!(read.G, 20);
        assert_eq!(read.N, 20);
    }
}

impl ColBaseComp {

    pub fn new (pos: usize) -> ColBaseComp {
        ColBaseComp {pos: pos, A: 0, T: 0, G: 0, C: 0, N: 0}
    }

    pub fn extract (&mut self, s: &u8) {
        match s {
            b'A' => self.A += 1,
            b'T' => self.T += 1,
            b'G' => self.G += 1,
            b'C' => self.C += 1,
            b'N' => self.N += 1,
            _ => panic!("Invalid character {:?} == {:?} found in read", *s, s.to_ascii_lowercase())
        }            
    }

    pub fn percentage (&mut self) {
        let sum = (self.A + self.G + self.C + self.T) as f64;
        self.A = ((self.A as f64 / sum)  * 100f64).round() as usize;
        self.G = ((self.G as f64 / sum)  * 100f64).round() as usize;
        self.T = ((self.T as f64 / sum)  * 100f64).round() as usize;
        self.C = ((self.C as f64 / sum)  * 100f64).round() as usize;
    }
}

#[derive(Serialize, Deserialize)]
pub struct BaseComp {
    lib: Vec<ColBaseComp>,
    len: usize,
}

impl BaseComp {
    pub fn init (len: usize) -> BaseComp {
        let mut base_comp = BaseComp { lib: Vec::with_capacity(len), len: 0 };
        for i in 0..len {
            base_comp.lib[i] = ColBaseComp::new(len);
        }

        base_comp
    }

    pub fn extract (&mut self, s: &str) {
        for c in s.as_bytes().iter().zip(0..s.len()) {
            self.lib[c.1].extract(c.0);
        }
    }

    pub fn jsonify (&self) -> String {
        serde_json::to_string(&self).expect("Error encountered JSONifying data")
    }
}