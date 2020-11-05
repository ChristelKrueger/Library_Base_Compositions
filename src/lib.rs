pub mod fastq_io {
    use std::path::PathBuf;
    use std::io::{self, BufRead, BufWriter, Write, BufReader, ErrorKind};
    use std::fs::OpenOptions;
    use std::process;

    /// Module to provide functionality for easy testing
    pub mod test_utils {
        use std::io::Cursor;
        /// Returns reader which implements Read trait. Reader yields data passed to parameter _s_
        pub fn return_reader<'a> (s: &'a[u8]) -> Cursor<&[u8]> {
           Cursor::new(s)
        }

        /** Returns a writer which implements Write trait.
        writer.get_ref()[0..].to_vec() can be used to get the data written to the writer.
        */
        pub fn return_writer() -> Cursor<Vec<u8>> {
            Cursor::new(Vec::<u8>::new())
        }
    }

    // impl read and write logic for FASTQRead struct
    use crate::sample::FASTQRead;
    impl FASTQRead {
        /**
         * Reads a complete FASTQ statement (composed of 4 lines) into the FASTQRead struct from given reader
         */
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


        /**
         * Prints the contents of FASTQRead. Pass false for a flag if it shouldn't be printed
         */
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

        #[test]
        fn test_percentage () {
            assert_eq!(run_with_content(b"AA\nTA\nGA\nCA\nNA\n"),
                std::str::from_utf8(br#"2
[{"pos":1,"A":20,"T":20,"G":20,"C":20,"N":20},{"pos":2,"A":100,"T":0,"G":0,"C":0,"N":0}]"#).unwrap())
        }

        use crate::fastq_io;
   
        fn run_with_content (content: &[u8]) -> String {
            let reader = fastq_io::test_utils::return_reader(content);
            let mut writer = fastq_io::test_utils::return_writer();

            run(reader, &mut writer);
            let stringified = writer.get_ref()[0..].to_vec();

            std::str::from_utf8(&stringified).unwrap().to_string()
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
    #[derive(Serialize, Deserialize, Clone, Debug, Copy, PartialEq)]
    #[allow(non_snake_case)]
    pub struct Read {
        pub pos: usize,
        pub A: usize,
        pub T: usize,
        pub G: usize,
        pub C: usize,
        pub N: usize,
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

        // Useful for making calculations on all 5 fields
        /// Returns all 5 bases as an array
        pub fn as_arr (&self) -> [usize; 5] {
            [self.A, self.T, self.G, self.C, self.N]
        }

        pub fn from_iter<I> (mut iter: I, pos: usize) -> Read
        where
            I: Iterator<Item = usize>,
        {
            let mut arr: [usize; 5] = [0; 5];
            for a in arr.iter_mut() {
                match iter.next() {
                    Some(n) => *a = n,
                    None => panic!("Improper iterator recieved, iterator is only supposed to return 5 elements for A, T, G, C, N respectively")
                }
            }

            Read {A: arr[0], T: arr[1], G: arr[2], C: arr[3], N: arr[4], pos: pos}
        }

        pub fn calc_percentage(&mut self) {
            let sum = self.as_arr().iter().fold(0.0, |acc, curr| acc + (*curr as f32));
            *self = Read::from_iter(self.as_arr().iter().map(|&x| ((x as f32 / sum) * 100.0).round() as usize), self.pos);
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
        
        //Convert raw numbers into %age and convert pos to start from 1
        for c in comp.iter_mut() {
            c.calc_percentage();
            c.pos += 1;
        }
        write!(writer, "{}\n", serde_json::to_string(&comp.len()).unwrap()).unwrap();
        write!(writer, "{}", serde_json::to_string(&comp).unwrap()).unwrap();
    }
}

pub mod plot_comp {
    use crate::base_extraction::Read;
    use serde_json;
    use structopt::StructOpt;
    use std::path::PathBuf;
    use plotters::prelude::*;

    #[derive(Debug, StructOpt)]
    #[structopt(name = "Plot base composition", about = "Plots base composition of given JSON file")]
    pub struct Cli {
        #[structopt(short = "I", long = "stdin")]
        /// Toggle stdin input
        stdin: bool,

        #[structopt(parse(from_os_str), required_unless("stdin"))]
        pub input: Option<PathBuf>,
        /// Other files of library for comparison
        #[structopt(short = "l", parse(from_os_str))]
        pub libs: Option<Vec<PathBuf>>,
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn test_calc_mean() {
            assert_eq!(calc_mean(
            &vec![
                vec![
                    Read {pos: 1, A: 7, T: 8, G: 55, C: 27, N: 2},
                ],
                vec![
                    Read {pos: 1, A: 7, T: 8, G: 53, C: 30, N: 2},
                ]
            ], 0), 
            Read {pos: 1, A: 7, T: 8, G: 54, C: 28, N: 2}
            )
        }

        #[test]
        fn test_calc_sd () {
            assert_eq!(calc_sd(&vec![
                vec![
                    Read {pos: 0, A: 25, T: 0, G: 75, C: 0, N: 10},
                ],
                vec![
                    Read {pos: 0, A: 75, T: 100, G: 100, C: 0, N: 10},
            ]], Read {pos: 0, A: 50, T: 50, G: 87, C: 0, N: 10},
             0), //This is mean of above values 
            Read {pos: 0, A: 25, T: 50, G: 12, C: 0, N: 0}
            )
        }
    }

    use std::io::BufRead;
    pub fn read_comp_file <R>(reader: &mut R) -> (usize, Vec<Read>)
    where R: BufRead {
        let mut s = String::with_capacity(500);

        //Read x_len from first line
        reader.read_line(&mut s).expect("Error reading line");
        s = s.trim().to_string();
        let x_len: usize = s.parse::<usize>().expect("Error reading max len, rerun extract-comp again");
        s.clear();

        //Read next line for JSON data
        reader.read_line(&mut s).expect("Error reading line");
        let mut comp: Vec<Read> = serde_json::from_str(&s).expect("Error converting JSON to data");
        comp.sort_by(|a, b| a.pos.cmp(&b.pos));

        (x_len, comp)
    }
    
    fn calc_mean (libs: &Vec<Vec<Read>>, pos: usize) -> Read {
        Read::from_iter(libs.iter()
        .map(move |lib| lib[pos]) //Get reads at pos of lib
        .fold(Read::new(pos), // Initial value: new read
                              // Adds values from reads (produced by map |lib| libs[pos]) to accumulated Read
                              // Then divides values of final accumulated read by x_len
            |acc, curr| Read::from_iter(
                acc.as_arr().iter() // Convert accumulated Read to its elements' iterator
                .zip(curr.as_arr().iter()) // Zip with current Read's elements' iterator                            
                .map(|x| x.0 + x.1) // Add their corresponding elements
            , pos)) 
        .as_arr().iter()
        .map(|x| x / libs.len()), pos + 1)
    }

    fn calc_sd (libs: &Vec<Vec<Read>>, mean: Read, pos: usize) -> Read {
        Read::from_iter(libs.iter()
            .map(move |lib| lib[pos]) //Get reads at pos of lib
            // Get differences from mean
            .map(|read| Read::from_iter(
                read.as_arr().iter() // Convert read to iterator over its elements
                .zip(mean.as_arr().iter()) // Zip it with iterator over mean's elements
                .map(|x| (*x.0 as isize - *x.1 as isize).abs() as usize) // Get difference b/w mean and element
                .map(|x| x * x), pos + 1)) // Square difference
            // Read::from_iter then converts it back into a Read
            
            // Get average of differences (squared)
            .fold(Read::new(pos + 1),
                |acc, curr| Read::from_iter(
                    acc.as_arr().iter() //Convert accumulated Read to its elements' iterator
                    .zip(curr.as_arr().iter()) //Zip with current Read's elements' iterator
                    .map(|x| x.0 + x.1) // Add the two numbers
                , pos)) //Calculates sum
            .as_arr().iter()
            .map(|x| x / libs.len()) // Get average
            .map(|x| (x as f64).sqrt() as usize) // Get square root of average
        //Read::from_iter turns the iterators of elements back into a read

        , pos)
    }

    use crate::fastq_io::get_reader;
    pub fn run <R>(mut reader: R, libs: Option<Vec<PathBuf>>) -> Result<(), Box<dyn std::error::Error>>
    where
        R: BufRead,
    {
        let (x_len, comp) = read_comp_file(&mut reader);

        //Set up plotting logic
        let root = BitMapBackend::new("out.png", (1024, 768)).into_drawing_area();
        root.fill(&WHITE)?;
        let mut chart = ChartBuilder::on(&root)
            // Set the caption of the chart
            .caption("Percentage", ("sans-serif", 40).into_font())
            // Set the size of the label region
            .x_label_area_size(30)
            .y_label_area_size(30)
            .margin(30)
            // Finally attach a coordinate on the drawing area and make a chart context
            .build_cartesian_2d(1usize..x_len, 0usize..100usize)?;

        chart.configure_mesh()
            .x_desc("Base number")
            .y_desc("Occurences")
            .axis_desc_style(("sans-serif", 15))
            .draw()?;

        //Draw bases
        //To customize colors, use ```&RGBColor(0, 0, 0)``` inplace of ```&MAGENTA```, etc.
        for i in [
            (&MAGENTA, "Base A"),
            (&BLUE, "Base T"),
            (&GREEN, "Base G"),
            (&CYAN, "Base C"),
            (&RED, "Unknown Base"),
        ].iter().zip([
            |r: &Read| r.A,
            |r: &Read| r.T,
            |r: &Read| r.G,
            |r: &Read| r.C,
            |r: &Read| r.N
        ].iter()) {
            chart
                .draw_series(LineSeries::new(
                    comp.iter().map(|r| (r.pos, i.1(r))),
                    i.0.0,
                ))?
                .label(i.0.1)
                .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], i.0.0));
        }

        if let Some(libs) = libs {
            let libs: Vec<Vec<Read>> = libs.iter()
                .map(|x| get_reader(&Some(x.to_path_buf())))
                .map(|mut x| read_comp_file(&mut x).1)
                .collect();

            let mut mean: Vec<Read> = Vec::with_capacity(x_len);
            let mut sd: Vec<Read> = Vec::with_capacity(x_len);

            println!("x_len: {}", x_len);
            for pos in 0..x_len {
                mean.push(calc_mean(&libs, pos));
                sd.push(calc_sd(&libs, mean[pos], pos));
                println!("Data was: {:?}, {:?}", libs[0][pos], libs[1][pos]);
                println!("mean is: {:?}. sd is: {:?}", mean[pos], sd[pos]);
            }

            for i in [
                &MAGENTA,
                &BLUE,
                &GREEN,
                &CYAN,
                &RED,
            ].iter().zip([
                |r: &Read| r.A,
                |r: &Read| r.T,
                |r: &Read| r.G,
                |r: &Read| r.C,
                |r: &Read| r.N
            ].iter()) {
                chart
                    .draw_series(std::iter::once(Polygon::new(
                        (0..x_len)
                            .map(|x| (x + 1, (i.1)(&mean[x]) as isize - (i.1)(&sd[x]) as isize))
                            .map(|x| (x.0, if x.1 < 0 {0usize} else {x.1 as usize})) //Round off may cause -ve value, so replace them with 0
                            .chain(
                            (0..x_len).map(|x| (x + 1, (i.1)(&mean[x]) + (i.1)(&sd[x]))))
                        .collect::<Vec<(usize, usize)>>(),
                        &i.0.mix(0.2)
                    )))?;
            }
        }

        chart
            .configure_series_labels()
            .background_style(&WHITE.mix(0.8))
            .border_style(&BLACK)
            .draw()?;

        Ok(())
    }
}
