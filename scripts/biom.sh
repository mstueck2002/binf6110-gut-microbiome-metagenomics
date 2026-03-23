#!/bin/bash
#SBATCH --job-name=kraken_biom
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=kraken_biom.out

#Set paths
REPORT_DIR="$SCRATCH/bracken_output"
OUT_FILE="$SCRATCH/bracken.biom"

#Run kraken-biom from source
$HOME/biom_env/bin/kraken-biom \
  "$REPORT_DIR"/*_bracken_species.report \
  -o "$OUT_FILE" \
  --fmt json
