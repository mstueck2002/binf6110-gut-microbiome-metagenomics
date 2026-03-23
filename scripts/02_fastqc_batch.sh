#!/bin/bash
#SBATCH --job-name=fastqc_array
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --array=1-147
#SBATCH --output=fastqc_%A_%a.out

#load necessary modules
module load fastqc

#set paths
FASTQ_DIR=$SCRATCH/fastq_files
OUT_DIR=$SCRATCH/fastqc_results

#generate output directory
mkdir -p $OUT_DIR

#specify files in array
FILE=$(ls $FASTQ_DIR/*.fastq | sed -n "${SLURM_ARRAY_TASK_ID}p")

#run fastqc using parallel array
fastqc -o $OUT_DIR $FILE
