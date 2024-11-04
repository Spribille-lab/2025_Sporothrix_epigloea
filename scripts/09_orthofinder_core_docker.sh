#!/bin/bash
#SBATCH --account=def-xxxxxxx
#SBATCH --time=3-0:00
#SBATCH --cpus-per-task=64
#SBATCH --job-name=ortho_msa
#SBATCH --output=orthofinder.script.core.Aug7.logs.out
#SBATCH --mem=249G


module load StdEnv/2020
module load apptainer/1.1.8

cd $SLURM_TMPDIR
cp -r $SCRATCH/2023_02_Sporothrix/09_orthology/orthofinder_input_core .
ls .
# orthofinder with iqtree installed
apptainer run --bind $(pwd):/input $SCRATCH/docker_builds/orthofinder_latest.sif orthofinder -t 64 -a 64 -M msa -A mafft -o /input/orthofinder_output_core_Aug7 -f /input/orthofinder_input_core
ls .
ls $SLURM_TMPDIR
cp -r $SLURM_TMPDIR/orthofinder_output_core_Aug7 $SCRATCH/2023_02_Sporothrix/09_orthology/.