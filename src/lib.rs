pub mod extract_comp;
pub mod plot_comp;

pub mod test_utils {
    use std::io::Cursor;
    /// Returns reader which implements Read trait.
    /// - `s`: Data which should be yielded by the reader upon read
    pub fn return_reader<'a> (s: &'a[u8]) -> Cursor<&[u8]> {
        
        Cursor::new(s)
    }

    /// Returns a writer which implements Write trait.
    /// `writer.get_ref()[0..].to_vec()` can be used to get the data written to the writer.
    pub fn return_writer() -> Cursor<Vec<u8>> {
        Cursor::new(Vec::<u8>::new())
    }

    pub fn get_writer_content(writer: Cursor<Vec<u8>>) -> String {
        std::str::from_utf8(&writer.get_ref()[0..]).unwrap().to_string()
    }
}

pub mod io_utils {
    use std::path::PathBuf;
    use std::fs::OpenOptions;
    use std::io::{self, BufReader, BufRead, Write};

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

    use std::io::ErrorKind;
    /// Will return writer to File if PathBuf can be opened, will panic if File unavailable
    /// And return writer to stdout if PathBuf not given
    pub fn get_writer(output: &Option<PathBuf>) -> Box<dyn Write> {
        match output {
            Some(file) => Box::new(OpenOptions::new().append(true).open(&file).unwrap_or_else(|error| {
                if error.kind() == ErrorKind::NotFound {
                    OpenOptions::new().create(true).write(true).open(&file).expect("Problem creating the file!")
                } else {
                    panic!("Problem opening the file: {:?}", error);
                }
            })),
            None => {
                Box::new(io::stdout())
            },
        }
    }
}

use serde::{Serialize, Deserialize};
use serde_json;

#[derive(Serialize, Deserialize, PartialEq, Debug)]
#[allow(non_snake_case)]
pub struct BaseCompCol {
    pub pos: usize,
    pub bases: BaseCompColBases,
}
#[derive(Serialize, Deserialize, PartialEq, Debug, Copy, Clone)]
#[allow(non_snake_case)]
pub struct BaseCompColBases {
    pub A: usize,
    pub T: usize,
    pub G: usize,
    pub C: usize,
    pub N: usize,
}

pub struct IterableBaseCompColBases {
    bases: BaseCompColBases,
    state: BaseCompColBasesIteratorState,
}

impl Iterator for IterableBaseCompColBases {
    type Item = usize;

    // VERY VERY IMP. If you change this, then also change the impl FromIterator<usize> for BaseCompColBases
    fn next(&mut self) -> Option<usize> {
        match self.state {
            BaseCompColBasesIteratorState::A => {self.state = BaseCompColBasesIteratorState::C; Some(self.bases.A)},
            BaseCompColBasesIteratorState::C => {self.state = BaseCompColBasesIteratorState::G; Some(self.bases.C)},
            BaseCompColBasesIteratorState::G => {self.state = BaseCompColBasesIteratorState::T; Some(self.bases.G)},
            BaseCompColBasesIteratorState::T => {self.state = BaseCompColBasesIteratorState::N; Some(self.bases.T)},
            BaseCompColBasesIteratorState::N => {self.state = BaseCompColBasesIteratorState::End; Some(self.bases.N)},
            _ => None
        }
    }
}

impl From<BaseCompColBases> for IterableBaseCompColBases {
    fn from(base_comp_col_bases: BaseCompColBases) -> IterableBaseCompColBases {
        IterableBaseCompColBases {
            bases: base_comp_col_bases,
            state: BaseCompColBasesIteratorState::A
        }
    }
}

#[derive(Clone, Copy)]
pub enum BaseCompColBasesIteratorState {
    A,
    C,
    G,
    T,
    N,
    End,
}

impl BaseCompColBases {
    pub fn iter(self) -> impl Iterator<Item = usize> {
        IntoIterator::into_iter(IterableBaseCompColBases::from(self))
    }

    pub fn new() -> BaseCompColBases {
        BaseCompColBases {A: 0, G: 0, T: 0, C: 0, N: 0}
    }

    pub fn percentage (&mut self) {
        let sum = self.iter().fold(0, |acc, curr| acc + curr) as f64;
        *self = self.iter().map(|base| ((base as f64 / sum)  * 100f64).round() as usize).collect();
    }
}

use std::iter::FromIterator;
impl FromIterator<usize> for BaseCompColBases {
    fn from_iter<I: IntoIterator<Item=usize>>(iter: I) -> Self {
        let mut c = IterableBaseCompColBases::from(BaseCompColBases {A: 0, G: 0, C:0, T:0, N:0});
        let mut iter = iter.into_iter();

        for ele in &mut
        // VERY VERY IMP. If you change this, then also change the impl Iterator for IterableBaseCompColBases
        // so the corresponding bases are collected into their corresponding spaces
            [&mut c.bases.A,
             &mut c.bases.C,
             &mut c.bases.G,
             &mut c.bases.T, 
             &mut c.bases.N].iter_mut() {
            **ele = match iter.next() {
                Some(n) => n,
                None => panic!("Improper iterator recieved"),
            }
        }

        c.bases
    }
}

#[cfg(test)]
mod col_base_comp_tests {
    use super::*;
    #[test]
    fn test_extract() {
        let mut read = BaseCompCol::new(0);
        read.extract(&b'A');
        assert_eq!(read.bases.A, 1);

        read.extract(&b'C');
        assert_eq!(read.bases.C, 1);

        read.extract(&b'T');
        assert_eq!(read.bases.T, 1);

        read.extract(&b'G');
        assert_eq!(read.bases.G, 1);

        read.extract(&b'N');
        assert_eq!(read.bases.N, 1);
    }
    #[test]
    fn test_percentage() {
        let mut read = BaseCompCol::new(0);
        for c in "ACTGN".as_bytes().iter() {
            read.extract(c);
        }
        read.bases.percentage();

        println!("{:?}", read);

        assert_eq!(read.bases.A, 20, "Testing A");
        assert_eq!(read.bases.C, 20, "Testing C");
        assert_eq!(read.bases.T, 20, "Testing T");
        assert_eq!(read.bases.G, 20, "Testing G");
        assert_eq!(read.bases.N, 20, "Testing N");
    }
}

impl BaseCompCol {
    pub fn new (pos: usize) -> BaseCompCol {
        BaseCompCol {pos: pos, bases: BaseCompColBases {A: 0, T: 0, G: 0, C: 0, N: 0}}
    }

    pub fn extract (&mut self, s: &u8) {
        match s {
            b'A' => self.bases.A += 1,
            b'T' => self.bases.T += 1,
            b'G' => self.bases.G += 1,
            b'C' => self.bases.C += 1,
            b'N' => self.bases.N += 1,
            _ => panic!("Invalid character {:?} == {:?} found in read", *s as char, s.to_ascii_lowercase())
        }            
    }

    
}

#[derive(Serialize, Deserialize, Debug)]
pub struct BaseComp {
    pub lib: Vec<BaseCompCol>,
    pub len: usize,
}

impl BaseComp {
    pub fn init (len: usize) -> BaseComp {
        let mut base_comp = BaseComp { lib: Vec::with_capacity(len), len};
        for i in 0..len {
            base_comp.lib.push(BaseCompCol::new(i + 1));
        }

        base_comp
    }

    pub fn extract (&mut self, s: &str) {
        for c in s.as_bytes().iter().enumerate() {
            self.lib[c.0].extract(c.1);
        }
    }

    pub fn jsonify (&self) -> String {
        serde_json::to_string(&self).expect("Error encountered JSONifying data")
    }
}