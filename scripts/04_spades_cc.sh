#!/bin/bash
#SBATCH --account=def-xxxxxxx
#SBATCH --time=3-0:0
#SBATCH --cpus-per-task=64
#SBATCH --mem=2000G
#SBATCH --job-name=SPAdes
#SBATCH --output=spades_2000G.script.logs.out


module load StdEnv/2020
module load spades/3.15.4

SAMPLE="$1" #ca202112
READS_1="$SAMPLE"_1.fastq.gz
READS_2="$SAMPLE"_2.fastq.gz

mkdir $SCRATCH/2023_01_Tmes_metagenomes/04_metagenome_assemblies/"$SAMPLE"
spades.py -t ${SLURM_CPUS_PER_TASK} -m 2000 --meta --pe1-1 $SCRATCH/2023_01_Tmes_metagenomes/03_CLEAN_READS/"$READS_1" --pe1-2 $SCRATCH/2023_01_Tmes_metagenomes/03_CLEAN_READS/"$READS_2" -k 21,31,51,71,81,101,127 -o $SCRATCH/2023_01_Tmes_metagenomes/04_metagenome_assemblies/"$SAMPLE"
