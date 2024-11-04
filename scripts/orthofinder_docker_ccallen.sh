#!/bin/bash
#SBATCH --account=def-xxxxxxxx
#SBATCH --time=3-0:00
#SBATCH --cpus-per-task=64
#SBATCH --job-name=ortho_msa
#SBATCH --output=orthofinder.script.Aug2.logs.out
#SBATCH --mem=249G

module load StdEnv/2020
module load apptainer/1.1.8

cd $SLURM_TMPDIR
cp -r $SCRATCH/2023_02_Sporothrix/09_orthology/orthofinder_input .
ls .
# orthofinder with iqtree installed
apptainer run --bind $(pwd):/input $SCRATCH/docker_builds/orthofinder_latest.sif orthofinder -t 64 -a 64 -M msa -A mafft -o /input/orthofinder_output_Aug3 -f /input/orthofinder_input
ls .
ls $SLURM_TMPDIR
cp -r $SLURM_TMPDIR/orthofinder_output_Aug3 $SCRATCH/2023_02_Sporothrix/09_orthology/.