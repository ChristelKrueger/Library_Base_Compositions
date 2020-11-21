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
            std::str::from_utf8(br#"{"lib":[{"pos":1,"A":20,"T":20,"G":20,"C":20,"N":20},{"pos":2,"A":100,"T":0,"G":0,"C":0,"N":0}],"len":2}"#).unwrap())
    }

    use fastq_io;

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

/// Can be used to send Vec<Read> to another program in JSON format
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct LibReads {
    pub lib: Vec<Read>,
    pub len: usize
}

impl Read {
    pub fn new (pos: usize) -> Read {
        Read {pos: pos, A: 0, T: 0, G: 0, C: 0, N: 0}
    }

    fn extract(&mut self, s: &u8) {
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

impl LibReads {
    /// Reads number of occurences per base into array of Reads, per column
    /// - `comp`: Array of Reads into which composition will be read. 
    /// Extra elements may be added if length of a line of FASTQ quality data is greater than length of given array.
    /// Occurences per base for first column will be read to first element of `comp`, and so on.
    /// - `reader`: Object implementing `std::io::BufRead` from which to read lines
    fn read <R>(&mut self, reader: &mut R) -> bool
    where
        R: BufRead
    {
        let mut s = String::with_capacity(self.lib.len());
        if reader.read_line(&mut s).expect("Error reading line. Make sure terminal supports UTF-8 input") == 0 {
            return false;
        }
        s = s.trim().to_string();

        // Initialize extra Reads to array if string length exceeds Read array size
        for i in self.lib.len()..s.len() {
            self.lib.push(Read::new(i))
        };
        for c in s.as_bytes().iter().zip(0..s.len()) {
            self.lib[c.1].extract(c.0);
        }
        true
    }

    fn update (&mut self) {
        self.len = self.lib.len();
    }
}

pub fn run<R, W>(mut reader: R, mut writer: W)
where
    R: BufRead,
    W: Write,
{
    // Each element represents a column of the seqs
    let mut comp: LibReads = LibReads{ lib: Vec::with_capacity(50), len: 0 }; //Won't have to resize array every time

    for c in comp.lib.iter_mut().zip(0..comp.len) {
        *c.0 = Read {pos: c.1, A: 0, T: 0, G: 0, C: 0, N: 0}; //Init all elements
    }

    //Read seqs into Read arr
    while LibReads::read(&mut comp, &mut reader) {}
    comp.update();
    
    //Convert raw numbers into %age and convert pos to start from 1
    for c in comp.lib.iter_mut() {
        c.calc_percentage();
        c.pos += 1;
    }

    comp.update();
    write!(writer, "{}", serde_json::to_string(&comp).unwrap()).unwrap();
}
