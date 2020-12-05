use fastq2comp::plot_comp::{run};
use fastq2comp::io_utils;
use std::path::PathBuf;
use plotters::prelude::BitMapBackend;

fn main() {
    run(
        io_utils::get_reader(&Some(PathBuf::from("sample_in.json"))),
        None,
        BitMapBackend::new(&PathBuf::from("sample_out.png"),
        (1280, 700))
    ).expect("Error drawing chart");
}
