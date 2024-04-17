#!/bin/bash
#SBATCH --account=def-tspribi
#SBATCH --time=3-0:00
#SBATCH --cpus-per-task=64
#SBATCH --job-name=ortho_msa
#SBATCH --output=orthofinder.script.Aug2.logs.out
#SBATCH --mem=249G
#SBATCH --mail-user=w6p9c9j6t9c6a2i6@spribillelabworkspace.slack.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

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