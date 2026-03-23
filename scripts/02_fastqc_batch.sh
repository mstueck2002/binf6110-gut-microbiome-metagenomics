#!/bin/bash
#SBATCH --job-name=fastqc_array
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --array=1-147
#SBATCH --output=fastqc_%A_%a.out

module load fastqc

FASTQ_DIR=$SCRATCH/fastq_files
OUT_DIR=$SCRATCH/fastqc_results

mkdir -p $OUT_DIR

FILE=$(ls $FASTQ_DIR/*.fastq | sed -n "${SLURM_ARRAY_TASK_ID}p")

fastqc -o $OUT_DIR $FILE
