use fastq2comp::plot_comp::{Cli, run};
use fastq2comp::io_utils;

use simple_logger::SimpleLogger;
use structopt::StructOpt;
use std::ffi::OsString;
use std::path::PathBuf;
use plotters::prelude::BitMapBackend;

fn main() {
    //Set up logger.
    SimpleLogger::new().init().unwrap();

    let args = Cli::from_args();

    // Output file logic
    let mut out_path: PathBuf = PathBuf::from(&args.input);
    let mut file_name = OsString::from(
        match out_path.file_stem() {
            Some(n) => n,
            None => panic!("Improper path provided"),
        }
    );
    file_name.push("_comp.png");

    out_path.set_file_name (file_name);

    run(io_utils::get_reader(&Some(args.input)), args.libs, BitMapBackend::new(&out_path, (1280, 700))).expect("Error drawing chart");
    eprintln!("Chart plotted at: {:#?}", &out_path);
}
