#!/bin/bash
# This script takes in n-number of fastq files and samples them, ready for plotting

for arg; do
  echo $arg
  cargo run --bin extract-fastq-comp 10000 $arg $arg"comp.json" 
done