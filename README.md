![Build](https://github.com/ChristelKrueger/Library_Base_Compositions/workflows/Rust/badge.svg)

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
*	[ ] Make nice front end/ability to upload own data
  * Babraham website
  * Online app?

# Library Base Compositions
Cambiohack project to make a QC tool to check for sequencing library base compositions

# Install Instructions:

### Setting up Sample SRR
* Install `Python 3` 

### Setting up composition extraction
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


### Setting up Metadata insertion in base composition file

* Install `jq` with `sudo apt install jq`.

### Final command:
```bash
# Substitute GDS_OUT_LOCATION with output file from the perl script.

cat data/$GDS_OUT_LOCATION | bash -x data/download-extract/download-extract.sh
```