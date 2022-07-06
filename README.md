[![Build & Test Rust](https://github.com/ChristelKrueger/Library_Base_Compositions/actions/workflows/test.yml/badge.svg)](https://github.com/ChristelKrueger/Library_Base_Compositions/actions/workflows/test.yml)

# Please note:
This repository is no longer active and its content has been integrated into https://github.com/DesmondWillowbrook/Librarian

# Problem Statement:
The base composition of sequencing reads depends on the library type (RNA, genomic, bisulfite, ChIP, etc.) and the species, and can often be characteristic for a particular sequencing application. For a while we’ve been thinking about a quality control tool that checks if a given base composition matches the expected base composition for the application. In other words, does my library look like it is supposed to? Some of the code of my last year’s hackathon project (Charades) could easily be adapted to put a given base composition into the wider context, but what’s missing is a collection of base compositions for a variety of sequencing libraries. The immediate task would be to think about how to best collect library base compositions and match them up with meta data about library type for a variety of published applications.

# Aims:
* Check base composition of your files match what you expect from that type of library
* Make the decision of the library selection on base composition of the library
  * Raise a red flag if there is no match before start of analysis 
# Tasks:
* [X] Compiling a set of commonly used sequencing libraries for bulk and single cell sequencing
  * Finding examples for these libraries (sequencing file plus metadata)
	* GEO accession number
		 * EBI ENA website?
* [X] Downloading 100000 reads (from FastQ file) for those libraries and storing them in a sensible way
*	[X] Extract base compositions per base
*	[X] Plot compositions
* [X] Sampling reads randomly instead of top _n_ reads
* [X] Make nice front end/ability to upload own data (Further development taking place in [Librarian](https://github.com/DesmondWillowbrook/Librarian-Server))
  * Babraham website
  * Online app?

# Library Base Compositions
Cambiohack project to make a QC tool to check for sequencing library base compositions

# Install Instructions:

### Setting up Sample SRR
* Install `Python 3` 

### Setting up composition extraction

* Install build dependencies, such as:
```
apt-get install cmake libfreetype6-dev libfontconfig1-dev pkg-config cmake build-essential
```
(Can possibly substitute `build-essential` with `clang`)

* Install Rust using `curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh`

* Make sure you have `Cargo` installed (can verify using `cargo --version`)

* Make sure latest version of Rust is installed. (run `rustup update`)

* Clone this repository:
   ```
   git clone https://github.com/ChristelKrueger/Library_Base_Compositions.git
   ```
   
* Then run:
   ```
   cargo build --release
   ```
   
* Your binaries will be in `target/release/`

### Final command:
```bash
# pwd should be root of the project, where this README is stored.

# Substitute GDS_OUT_LOCATION with output file from the perl script.
bash data/download-extract/download-extract.sh < data/$GDS_OUT_LOCATION

# In the case of the sample output, it would be:
bash data/download-extract/download-extract.sh < data/results_example_gds.txt
```
Output will be appended to `output.csv` file.
