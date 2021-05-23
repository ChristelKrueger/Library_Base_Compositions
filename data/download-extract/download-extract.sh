#!/bin/bash

# Should be run from project root.
# Test command: cat data/results_example_gds.txt | bash -x data/download-extract/download-extract.sh

# Made with reference to the following format:
# 1	Mus musculus	RNA-Seq	SRR12478073	ftp.sra.ebi.ac.uk/vol1/fastq/SRR124/073/SRR12478073/SRR12478073.fastq.gz	Transcriptome analysis of the TA muscles from WT and Dmd Exon 51 Knockout mice
# 2	Mus musculus	RNA-Seq	SRR12926516	ftp.sra.ebi.ac.uk/vol1/fastq/SRR129/016/SRR12926516/SRR12926516_1.fastq.gz	Prenatal exposure to environmental toxins induces sexually dimorphic transcriptome changes in neonatal offspring

# i.e with a \t (TAB) character between each bit of info

# So script doesn't keep going after error
set -euo pipefail

while IFS='$\n' read -r line; do
    srr_number=`echo "$line" | awk -F '\t' '{ print $4 }'`
    echo -e "$line\t$(python3 ./data/download-extract/sample_srr.py $srr_number 2 100 | ./target/release/extract_comp --stdin --stdout --tsv --trim 50 100)" >> ./data/download-extract/output.tsv
done
