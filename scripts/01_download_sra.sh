#!/bin/bash
#SBATCH --job-name=sra_download
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=sra_download.out

#load necessary modules
module load sra-toolkit

#prefetch SRR IDs
xargs -a srr_ids.txt -n 1 -P 4 prefetch

#find associated SRA and retrieve paired fastq files
find . \( -name "*.sra" -o -name "*.sralite" \) -print0 | \
xargs -0 -n 1 -P 2 fasterq-dump --split-files --threads 8
