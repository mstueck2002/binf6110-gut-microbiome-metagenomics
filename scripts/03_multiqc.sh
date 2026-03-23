#!/usr/bin/env bash
#Requires: conda env "multiqc_env" 

set -euo pipefail

#set paths
FASTQC=$SCRATCH/fastqc_results
OUT=$SCRATCH/multiqc_report

#generate output directory 
mkdir -p "$OUT"

#run multiqc on collective fastqc results
multiqc "$FASTQC" -o "$OUT"
