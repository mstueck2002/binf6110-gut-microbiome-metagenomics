#!/usr/bin/env bash
#Requires: conda env "multiqc_env" 

set -euo pipefail

FASTQC=$SCRATCH/fastqc_results
OUT=$SCRATCH/multiqc_report

mkdir -p "$OUT"

multiqc "$FASTQC" -o "$OUT"
