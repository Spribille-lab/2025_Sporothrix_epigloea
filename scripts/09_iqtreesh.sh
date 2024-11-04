#!/bin/bash
#SBATCH --account=def-xxxxxxx
#SBATCH --time=5-0:0
#SBATCH --cpus-per-task=29
#SBATCH --mem=200G
#SBATCH --job-name=iqtree
#SBATCH --output=iqtree.script.logs.out


module load StdEnv/2020
module load gcc/9.3.0
module load iq-tree/2.0.7

cd /home/ccallen/scratch/2023_02_Sporothrix/09_orthology/orthofinder_output_Aug3/Results_Aug03/Species_Tree/iqtree2
iqtree2 -s ../../MultipleSequenceAlignments/SpeciesTreeAlignment.fa -o Leptographium_lundbergii_CBS138716 -nt 29 -m MFP -B 1000 -pre sporothrix_placement