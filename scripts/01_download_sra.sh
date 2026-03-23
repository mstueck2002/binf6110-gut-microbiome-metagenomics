#!/bin/bash
#SBATCH --job-name=sra_download
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=sra_download.out

module load sra-toolkit

xargs -a srr_ids.txt -n 1 -P 4 prefetch

find . \( -name "*.sra" -o -name "*.sralite" \) -print0 | \
xargs -0 -n 1 -P 2 fasterq-dump --split-files --threads 8
