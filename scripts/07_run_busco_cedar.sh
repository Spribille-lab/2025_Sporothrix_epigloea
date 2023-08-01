#!/bin/bash
#SBATCH --account=def-tspribi
#SBATCH --time=2-0:0
#SBATCH --cpus-per-task=8
#SBATCH --job-name=array
#SBATCH --output=busco.run2.script.out
#SBATCH --mem=100G
#SBATCH --mail-user=w6p9c9j6t9c6a2i6@spribillelabworkspace.slack.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

bin_folder=$1

module load singularity

cd $SLURM_TMPDIR
cp $HOME/projects/def-tspribi/busco_db_2022/busco_downloads.tar.gz .
tar -x --use-compress-program="pigz -p 8 " -f busco_downloads.tar.gz

cd $SCRATCH/2023_01_Tmes_metagenomes/06_binning/"$bin_folder"

for file in $(ls *.fa)
do
echo the file is $file
genome=${file%.*}
echo the genome is $genome
singularity run --bind $(pwd):/data,$SLURM_TMPDIR/busco_downloads:/busco_downloads --pwd /data $SCRATCH/busco_v5.3.2_cv1.sif busco -i $file -o $genome -m genome --cpu 8 -f --download_path /busco_downloads --offline
#cleaning
find /scratch/ccallen/2023_01_Tmes_metagenomes/06_binning/"$bin_folder"/"$genome" -maxdepth 1 -mindepth 1 -type d -exec tar -czf {}.tar.gz {} --remove-files \;
done