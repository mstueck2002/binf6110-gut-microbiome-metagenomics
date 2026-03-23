#!/bin/bash
#SBATCH --job-name=bracken_array
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --array=1-74%8
#SBATCH --output=bracken_%A_%a.out

#Load necessary modules
module load bracken

#Set paths
REPORT_DIR="$SCRATCH/kraken_reports"
OUT_DIR="$SCRATCH/bracken_output"
DB_DIR="$SCRATCH/kraken_kb"
SAMPLE_LIST="$SCRATCH/sra_download/srr_ids.txt"

#Generate output directory
mkdir -p "$OUT_DIR"

#Get the sample ID for this array task
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

#Define input report file
REPORT_FILE="${REPORT_DIR}/${sample}.report"

#Run Bracken
bracken \
  -d "$DB_DIR" \
  -i "$REPORT_FILE" \
  -o "${OUT_DIR}/${sample}.bracken" \
  -l S \
  -r 150
