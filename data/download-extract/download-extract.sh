#!/bin/bash
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

    # Generates filter file for jq to use
    filter_file=$(mktemp)
    # $'\"' is escaped double quote ("). God, why does this need to be explaned? Because bash.
    echo {metadata: {serial_num:$'\"'$serial_num$'\"', species:$'\"'$species$'\"', lib_type:$'\"'$lib_type$'\"', srr_number:$'\"'$srr_number$'\"', title:$'\"'$title$'\"'}, data:.} > $filter_file

    wget -q -O - $ftp_url | gunzip -c | cargo run --bin extract_comp -- --stdin --stdout 100000 | \
    jq -f $filter_file > $srr_number"_comp.json"
done