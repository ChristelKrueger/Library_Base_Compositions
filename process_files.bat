rem This script takes in n-number of fastq files and samples them, ready for plotting

@echo off
setlocal enableDelayedExpansion
for %%x in (%*) do (
    echo %%x
    echo %%x_filtered.txt
    cargo run --bin sample-fastq 10000 %%x %%x_filtered.txt --skip-quals --skip-header --skip-mid
    cargo run --bin extract-comp %%x_filtered.txt > %%x_filtered_json.txt
)
endlocal
pause
    