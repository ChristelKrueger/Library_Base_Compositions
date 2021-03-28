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

    serial_num=`echo "$line" | awk -F '\t' '{ print $1 }'`
    species=`echo "$line" | awk -F '\t' '{ print $2 }'`
    lib_type=`echo "$line" | awk -F '\t' '{ print $3 }'`
    srr_number=`echo "$line" | awk -F '\t' '{ print $4 }'`
    ftp_url=`echo "$line" | awk -F '\t' '{ print $5 }'`
    title=`echo "$line" | awk -F '\t' '{ print $6 }'`
    METADATA="$serial_num,$species,$lib_type,$srr_number,$title"

    echo "$METADATA,$(python3 ./data/download-extract/sample_srr.py $srr_number 2 1000 | ./target/release/extract_comp --stdin --stdout --csv --trim 50 1000)" >> ./data/download-extract/output.csv
done