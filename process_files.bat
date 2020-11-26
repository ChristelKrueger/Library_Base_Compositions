rem This script takes in n-number of fastq files and samples them, ready for plotting

@echo off
setlocal enableDelayedExpansion
for %%x in (%*) do (
    echo %%x
    cargo run --bin extract-fastq-comp 10000 %%x %%x_comp.json
)
endlocal
pause
    