#!/bin/bash
# This script takes in n-number of fastq files and samples them, ready for plotting

for arg; do
  echo $arg
  echo $arg"_filtered.txt"
  cargo run --bin sample-fastq 10000 $arg $arg"_filtered.txt" --skip-quals --skip-header --skip-mid
  cargo run --bin extract-comp $arg"_filtered.txt" $arg"_filtered_json.txt"    
done