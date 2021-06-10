use fastq2comp::plot_comp::{run};
use fastq2comp::io_utils;
use std::path::PathBuf;
use plotters::prelude::BitMapBackend;

fn main() {
    run(
        io_utils::get_reader(&Some(PathBuf::from("examples/plot-comp/in.json")), false),
        None,
        BitMapBackend::new(&PathBuf::from("examples/plot-comp/out.png"),
        (1280, 700))
    ).expect("Error drawing chart");
}
