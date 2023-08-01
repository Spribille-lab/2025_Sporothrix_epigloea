#!/bin/bash
#SBATCH --account=def-tspribi
#SBATCH --time=2-0:00
#SBATCH --cpus-per-task=8
#SBATCH --job-name=ortho_msa
#SBATCH --output=orthofinder.script.logs.out
#SBATCH --mem=32G
#SBATCH --mail-user=w6p9c9j6t9c6a2i6@spribillelabworkspace.slack.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

module load StdEnv/2020
module load apptainer/1.1.8

cd $SLURM_TMPDIR
cp -r $SCRATCH/2023_02_Sporothrix/09_orthology/orthofinder_input .
ls .
# orthofinder with iqtree installed
apptainer run --bind $(pwd):/input $SCRATCH/docker_builds/orthofinder_latest.sif orthofinder -M msa -A mafft -T iqtree -o /input/orthofinder_output -f /input/orthofinder_input
ls .
ls $SLURM_TMPDIR
cp -r $SLURM_TMPDIR/orthofinder_output $SCRATCH/2023_02_Sporothrix/09_orthology/.

