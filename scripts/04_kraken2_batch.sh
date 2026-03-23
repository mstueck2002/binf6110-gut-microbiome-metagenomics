#!/bin/bash
#SBATCH --job-name=kraken2_ram
#SBATCH --cpus-per-task=16
#SBATCH --mem=0
#SBATCH --time=24:00:00
#SBATCH --array=1-6
#SBATCH --output=kraken_%A_%a.out

#load necessary modules
module load kraken2

#set database paths
DB_SRC=$SCRATCH/kraken_db
DB_RAM=/dev/shm/kraken_db

#set file paths
FASTQ=$SCRATCH/fastq_files
OUT=$SCRATCH/kraken_output
SAMPLE_LIST=$SCRATCH/sra_download/srr_ids.txt

#generate output directory
mkdir -p "$OUT"

#split work across array tasks
TOTAL_SAMPLES=$(wc -l < $SAMPLE_LIST)
TASKS=6

SAMPLES_PER_TASK=$(( (TOTAL_SAMPLES + TASKS - 1) / TASKS ))

START=$(( (SLURM_ARRAY_TASK_ID - 1) * SAMPLES_PER_TASK + 1 ))
END=$(( SLURM_ARRAY_TASK_ID * SAMPLES_PER_TASK ))

echo "Task $SLURM_ARRAY_TASK_ID processing samples $START to $END"

#copy the Kraken2 database to RAM (1x per job)
echo "Copying Kraken2 DB to RAM..."
mkdir -p $DB_RAM
rsync -a --info=progress2 $DB_SRC/ $DB_RAM/
echo "DB ready"

#process multiple samples per array
#define if/else statements for paired and the single-end read (73 paired, 1 single end)
#run kraken2 in parallel array
sed -n "${START},${END}p" $SAMPLE_LIST | while read sample
do
    echo "Processing sample: $sample"

    if [ -f "${FASTQ}/${sample}_2.fastq" ]; then
        kraken2 \
            --db "$DB_RAM" \
            --threads $SLURM_CPUS_PER_TASK \
	    --confidence 0.15 \
            --quick \
            --paired "${FASTQ}/${sample}_1.fastq" "${FASTQ}/${sample}_2.fastq" \
            --report "${OUT}/${sample}.report" \
            --output "${OUT}/${sample}.kraken"
    else
        kraken2 \
            --db "$DB_RAM" \
            --threads $SLURM_CPUS_PER_TASK \
	    --confidence 0.15 \
            --quick \
            "${FASTQ}/${sample}_1.fastq" \
            --report "${OUT}/${sample}.report" \
            --output "${OUT}/${sample}.kraken"
    fi

done

#Remove the kraken2 database once finished
rm -rf $DB_RAM
echo "Finished task $SLURM_ARRAY_TASK_ID"
