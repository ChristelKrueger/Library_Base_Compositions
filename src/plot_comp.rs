use crate::{
    BaseCompCol as Read,
    BaseComp as LibReads,
    BaseCompColBases as Bases,
    io_utils
};

use serde_json;
use structopt::StructOpt;
use std::path::PathBuf;
use plotters::prelude::*;


#[derive(Debug, StructOpt)]
#[structopt(name = "Plot base composition", about = "Plots base composition of given JSON file")]
pub struct Cli {
    #[structopt(parse(from_os_str))]
    pub input: PathBuf,
    /// Other files of library for comparison
    #[structopt(short = "l", parse(from_os_str))]
    pub libs: Option<Vec<PathBuf>>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::data_transforms::*;

    #[test]
    fn test_calc_mean() {
        assert_eq!(calc_mean(
        &vec![
            vec![
                Read {pos: 1, bases: Bases {A: 7, T: 8, G: 55, C: 27, N: 2}},
            ],
            vec![
                Read {pos: 1, bases: Bases {A: 7, T: 8, G: 53, C: 30, N: 2}},
            ]
        ], 0), 
        Bases {A: 7, T: 8, G: 54, C: 28, N: 2}
        )
    }

    #[test]
    fn test_calc_sd () {
        assert_eq!(calc_sd(&vec![
            vec![
                Read {pos: 1, bases: Bases {A: 25, T: 0, G: 75, C: 0, N: 10}},
            ],
            vec![
                Read {pos: 1, bases: Bases {A: 75, T: 100, G: 100, C: 0, N: 10}},
        ]], Bases {A: 50, T: 50, G: 87, C: 0, N: 10},
            0), //This is mean of above values 
        Bases {A: 35, T: 71, G: 18, C: 0, N: 0}
        )
    }
}

use std::io::BufRead;
pub fn read_comp_file <R>(reader: &mut R) -> (usize, Vec<Read>)
where R: BufRead {
    let mut s = String::with_capacity(500);

    //Read next line for JSON data
    reader.read_line(&mut s).expect("Error reading line");
    let mut comp: LibReads = serde_json::from_str(&s).expect("Error converting JSON to data");
    comp.lib.sort_by(|a, b| a.pos.cmp(&b.pos));

    (comp.len, comp.lib)
}

mod data_transforms {
    use super::*;

    pub fn calc_mean (libs: &Vec<Vec<Read>>, pos: usize) -> Bases {
        libs.iter()
        .map(move |lib| lib[pos].bases)
        .fold( Bases::new(),
            |acc, curr|
                acc.iter()
                .zip(curr.iter())
                .map(|base| base.0 + base.1)
                .collect()
        ) //Calculates sum of Bases
        .iter()
        .map(|x| x / libs.len()) //Divides each base with number of bases
        .collect()
    }

    pub fn calc_sd (libs: &Vec<Vec<Read>>, mean: Bases, pos: usize) -> Bases {
        libs.iter()
        .map(move |lib| lib[pos].bases)
        //Get differences from mean
        .map(|read|
            // Get difference b/w mean and element
            read.iter() // Convert read to iterator over its elements
            .zip(mean.iter()) // Zip it with iterator over mean's elements
            .map(|x| (x.0 as isize - x.1 as isize).abs() as usize) 

            // Square difference
            .map(|x| x * x) 
            .collect()
        )
        // Calculates sum of squared values in Bases
        .fold(Bases::new(),
            |acc, curr: Bases|
            acc.iter()
            .zip(curr.iter())
            .map(|base| base.0 + base.1)
            .collect()
        )
        .iter()
        // Divides squared values with length of libs to get average (i.e. variance)
        .map(|x| x / (if libs.len() > 1 {libs.len() - 1} else {1})) // Get average (subtract by 1 as this is sample data)
        .map(|x| (x as f64).sqrt().round() as usize) // Get square root of average
        .collect()
    }
}

use io_utils::get_reader;
use data_transforms::*;
use plotters_backend::DrawingBackend;

pub fn run <R, B>(mut reader: R, libs: Option<Vec<PathBuf>>, backend: B) -> Result<(), Box<dyn std::error::Error>>
where
    R: BufRead,
    B: DrawingBackend,
{
    let (x_len, comp) = read_comp_file(&mut reader);

    //Set up plotting logic
    let root = backend.into_drawing_area();
    root.fill(&WHITE).expect("Error drawing chart.");
    let mut chart = ChartBuilder::on(&root)
        // Set the caption of the chart
        .caption("Percentage", ("sans-serif", 40).into_font())
        // Set the size of the label region
        .x_label_area_size(30)
        .y_label_area_size(30)
        .margin(30)
        // Finally attach a coordinate on the drawing area and make a chart context
        .build_cartesian_2d(1usize..x_len, 0usize..100usize).expect("Error drawing chart.");

    chart.configure_mesh()
        .x_desc("Base number")
        .y_desc("Occurences")
        .axis_desc_style(("sans-serif", 15))
        .draw().expect("Error drawing chart.");

    //Draw bases
    //To customize colors, use ```&RGBColor(0, 0, 0)``` inplace of ```&MAGENTA```, etc.
    for i in [
        (&MAGENTA, "Base A"),
        (&BLUE, "Base T"),
        (&GREEN, "Base G"),
        (&CYAN, "Base C"),
        (&RED, "Unknown Base"),
    ].iter().zip([
        |r: &Read| r.bases.A,
        |r: &Read| r.bases.T,
        |r: &Read| r.bases.G,
        |r: &Read| r.bases.C,
        |r: &Read| r.bases.N
    ].iter()) {
        chart
            .draw_series(LineSeries::new(
                comp.iter().map(|r| (r.pos, i.1(r))),
                i.0.0,
            )).expect("Error drawing chart.")
            .label(i.0.1)
            .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], i.0.0));
    }

    if let Some(libs) = libs {
        let libs: Vec<Vec<Read>> = libs.iter()
            .map(|x| get_reader(&Some(x.to_path_buf())))
            .map(|mut x| read_comp_file(&mut x).1)
            .collect();

        for i in [
            &MAGENTA,
            &BLUE,
            &GREEN,
            &CYAN,
            &RED,
        ].iter().zip([
            |base: &Bases| base.A,
            |base: &Bases| base.T,
            |base: &Bases| base.G,
            |base: &Bases| base.C,
            |base: &Bases| base.N
        ].iter()) {
            chart
                .draw_series(
                    (0..x_len)
                        .filter_map(|pos| {
                            let mean = calc_mean(&libs, pos);
                            let sd = i.1(&calc_sd(&libs, mean, pos));
                            let mean = i.1(&mean);

                            // Logic to prevent overflows
                            // Negative results may be produced due to rounding-off
                            // So convert -ve results to 0
                            let lower = mean as isize - sd as isize;
                            let lower: usize = if lower < 0 {0} else {lower as usize};
                            let upper = mean + sd;

                            if upper != lower {
                                Some(ErrorBar::new_vertical(
                                    pos + 1,
                                    lower,
                                    mean,
                                    upper,
                                    (&i.0.mix(0.2)).filled(),
                                    5
                                ))
                            } else {
                                None
                            }
                        })                            
                ).expect("Error drawing chart.");
        }
    }

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw().expect("Error drawing chart.");

    Ok(())
}
