#!/bin/bash
#SBATCH --job-name=bracken_combine
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=bracken_combine.out

#Load necessary modules
module load bracken

#Set paths
OUTPUT_DIR="$SCRATCH/bracken_output"

#Run combine_bracken_outputs.py
python $EBROOTBRACKEN/analysis_scripts/combine_bracken_outputs.py \
  --files "$OUTPUT_DIR"/*.bracken \
  -o $SCRATCH/species_abundance_matrix.txt
