pub mod fastq_io {
    use std::path::PathBuf;
    use std::io::{self, BufRead, BufWriter, Write, BufReader, ErrorKind};
    use std::fs::OpenOptions;
    use std::process;

    pub mod test_utils {
        use std::io::Cursor;
        pub fn return_reader<'a> (s: &'a[u8]) -> Cursor<&[u8]> {
           Cursor::new(s)
        }

        pub fn return_writer() -> Cursor<Vec<u8>> {
            Cursor::new(Vec::<u8>::new())
        }
    }

    // impl read and write logic for FASTQRead struct
    use crate::sample::FASTQRead;
    impl FASTQRead {
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


        pub fn write<W> (&self, writer: &mut W, header: bool, seq: bool, mid: bool, quals: bool)
        where
            W: Write,
        {
            if header {write!(writer, "{}\n", self.header).expect("Error writing to file/stdout.");}
            if seq {write!(writer, "{}\n", self.seq).expect("Error writing to file/stdout.");}
            if mid {write!(writer, "{}\n", self.mid).expect("Error writing to file/stdout.");}
            if quals {write!(writer, "{}\n", self.quals).expect("Error writing to file/stdout.");}
        }
    }

    use crate::base_extraction::Read;
    impl Read {
        pub fn read <R>(comp: &mut Vec<Read>, reader: &mut R) -> bool
        where
            R: BufRead
        {
            let mut s = String::with_capacity(comp.len());
            if reader.read_line(&mut s).expect("Error reading line. Make sure terminal supports UTF-8 input") == 0 {
                return false;
            }
            s = s.trim().to_string();

            for i in comp.len()..s.len() {
                comp.push(Read::new(i))
            };
            for c in s.as_bytes().iter().zip(0..s.len()) {
                comp[c.1].extract(c.0);
            }
            true
        }

        pub fn read_json <R> (reader: &mut R) -> Option<Read>
        where 
            R: BufRead
        {
            let mut s = String::new();
            if reader.read_line(&mut s).expect("Error reading line. Make sure terminal supports UTF-8 input") == 0 {
                return None;
            }
            Some(serde_json::from_str(&s).expect("Error converting from JSON to data"))
        }
    }

    // Reader is a wrapper over BufRead
    // Takes in a PathBuf and open it or if no PathBuf is provided, opens up stdin
    // And provides an interface over the actual reading.
 
    pub fn get_reader(input: &Option<PathBuf>) -> Box<dyn BufRead> {
        match input {
            Some(file) => Box::new(BufReader::new(OpenOptions::new().read(true).open(&file).expect("Input file does not exist"))),
            None => {
                let stdin = io::stdin();
                Box::new(BufReader::new(stdin))
            },
        } 
    }

    pub fn get_writer(output: &Option<PathBuf>) -> BufWriter<Box<dyn Write>> {
        match output {
            Some(file) => BufWriter::new(Box::new(OpenOptions::new().append(true).open(&file).unwrap_or_else(|error| {
                if error.kind() == ErrorKind::NotFound {
                    OpenOptions::new().create(true).write(true).open(&file).expect("Problem creating the file!")
                } else {
                    panic!("Problem opening the file: {:?}", error);
                }
            }))),
            None => {
                BufWriter::new(Box::new(io::stdout()))
            },
        }
    }
}

pub mod sample {
    #[cfg(tests)]
    mod tests {
        use super::*;
        use super::test_utils::{return_reader, return_writer};

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

        use super::*;
        use super::test_utils::return_reader;

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

    pub fn run<R, W>(mut reader: R, mut writer: W, args: Cli)
    where
        R: BufRead,
        W: Write,
    {
        info!("Arguments recieved: {:#?}", args);

        //This is done as we want default value to be usize::MAX, which is not supported by StructOpt
        //This is a workaround
        let n_content = match args.n_content {Some(t) => t, None => usize::MAX};

        let mut read = FASTQRead::new(0);

        let mut valid_seqs: u64 = 0;

        // Set up write buffer to speed up write performance.
        // As we won't have to flush everytime we write.
        // TODO: There is an issue with BufWriter, it is not flushing automatically,
        // and flushing at the end does not flush the data in buffer inserted during the loop.
        // So as a stopgap we are flushing after printing a single read, which defeats the purpose of
        // using BufWriter

        loop {
            read.read(&mut reader);
            if let Some(n) = args.trimmed_length {
                read.trim(n);
            }

            // Check for numbers in reads
            if read.check_colorspace() {
                info!("Found numbers in reads - this is probably colorspace");
                break;
            }

            // Count the N's
            if read.count_n() > n_content {
                debug!("N count too high - skipping");
                continue;
            }

            if read.get_average_quality() < args.min_phred_score {
                debug!("Quality too low - skipping");
                continue;
            }

            read.write(&mut writer, !args.header, !args.seq, !args.mid, !args.quals);
            writer.flush().expect("Error flushing stream");

            valid_seqs += 1;

            if valid_seqs >= args.target_read_count {
                info!("Found enough sequences\n");
                break;
            }
        }
    }
}

pub mod base_extraction {
    #[cfg(test)]
    mod tests {
        use super::*;
        #[test]
        fn test_extract() {
            let mut read = Read::new(0);
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

        use crate::fastq_io;
        #[test]
        fn test_base_extraction() {
            let reader = fastq_io::test_utils::return_reader(b"AT\nGC\nNN\n");
            let mut writer = fastq_io::test_utils::return_writer();

            run(reader, &mut writer);
            let stringified = writer.get_ref()[0..].to_vec();
            assert_eq!(std::str::from_utf8(&stringified).unwrap().to_string(),
                std::str::from_utf8(
                br#"{"pos":0,"A":1,"T":0,"G":1,"C":0,"N":1}
{"pos":1,"A":0,"T":1,"G":0,"C":1,"N":1}
"#).unwrap())
        }
    }

    use structopt::StructOpt;
    use std::path::PathBuf;
    use serde::{Serialize, Deserialize};
    use serde_json;
    use std::io::{BufRead, Write};

    #[derive(Debug, StructOpt)]
    #[structopt(name = "extract fastq file", about = "Extracts base composition of given FASTQ file")]
    pub struct Cli {
        #[structopt(short = "I", long = "stdin")]
        /// Toggle stdin input
        stdin: bool,

        /// Input FASTQ file name
        #[structopt(parse(from_os_str), required_unless("stdin"))]
        pub input: Option<PathBuf>,
    }

    // Deriving Clone to allow easy initialization of Vec<Read>
    // pos included so multithreading will be easy since recieving program will be able to make sense of these
    // even if compositions arive out of order
    #[derive(Serialize, Deserialize, Clone)]
    #[allow(non_snake_case)]
    pub struct Read {
        pos: usize,
        A: u32,
        T: u32,
        G: u32,
        C: u32,
        N: u32,
    }

    impl Read {
        pub fn new (pos: usize) -> Read {
            Read {pos: pos, A: 0, T: 0, G: 0, C: 0, N: 0}
        }

        pub fn extract(&mut self, s: &u8) {
            match s {
                b'A' => self.A += 1,
                b'T' => self.T += 1,
                b'G' => self.G += 1,
                b'C' => self.C += 1,
                b'N' => self.N += 1,
                _ => panic!("Invalid character {:?} == {:?} found in read", *s, s.to_ascii_lowercase())
            }            
        }
    }
    
    pub fn run<R, W>(mut reader: R, mut writer: W)
    where
        R: BufRead,
        W: Write,
    {
        // Each element represents a column of the seqs
        let mut comp: Vec<Read> = Vec::new();
 
        let len = comp.len();
        for c in comp.iter_mut().zip(0..len) {
            *c.0 = Read {pos: c.1, A: 0, T: 0, G: 0, C: 0, N: 0};
        }

        //Read seqs into Read arr
        while Read::read(&mut comp, &mut reader){}
        write!(writer, "{}", serde_json::to_string(&comp).unwrap()).unwrap();
    }
}

pub mod draw {
    use crate::base_extraction::Read;
    use crate::fastq_io::get_reader;
    use serde_json::{self, Deserializer};
    use structopt::StructOpt;
    use std::path::PathBuf;

    #[derive(Debug, StructOpt)]
    #[structopt(name = "Plot base composition", about = "Plots base composition of given JSON file")]
    pub struct Cli {
        #[structopt(short = "I", long = "stdin")]
        /// Toggle stdin input
        stdin: bool,

        /// Input FASTQ file name
        #[structopt(parse(from_os_str), required_unless("stdin"))]
        pub input: Option<PathBuf>,
    }

    fn run(args: Cli) {
        let mut reader = get_reader(&args.input);
        let stream = Deserializer::from_reader(reader).into_iter::<Read>();
        
        for read in stream.into_iter() {
            // Plot a point on various lines for A, G, T, C
        }
    }
}