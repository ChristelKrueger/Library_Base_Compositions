use std::path::PathBuf;
use std::io::{self, BufRead, BufWriter, Write, BufReader, ErrorKind};
use std::fs::OpenOptions;

/// Module to provide functionality for easy testing
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

/// Will return writer to File if PathBuf can be opened, will panic if File unavailable
/// And return writer to stdout if PathBuf not given
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