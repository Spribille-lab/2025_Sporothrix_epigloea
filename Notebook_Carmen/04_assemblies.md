# assemble each Sporothrix genome using metaSPades 3.15.4 (cc)
# 2023.05.18


/home/ccallen/scratch/2023_02_Sporothrix/scripts/spades_cc.sh
```
#!/bin/bash
#SBATCH --account=def-tspribi
#SBATCH --time=3-0:0
#SBATCH --cpus-per-task=64
#SBATCH --mem=2000G
#SBATCH --job-name=SPAdes
#SBATCH --output=spades_2000G.script.logs.out
#SBATCH --mail-user=w6p9c9j6t9c6a2i6@spribillelabworkspace.slack.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

module load StdEnv/2020
module load spades/3.15.4

SAMPLE="$1"
READS_1="$SAMPLE"_1.fastq.gz
READS_2="$SAMPLE"_2.fastq.gz

mkdir $SCRATCH/2023_02_Sporothrix/04_assemblies/"$SAMPLE"
spades.py -t ${SLURM_CPUS_PER_TASK} -m 2000 --meta --pe1-1 $SCRATCH/2023_02_Sporothrix/03_CLEAN_READS/"$READS_1" --pe1-2 $SCRATCH/2023_02_Sporothrix/03_CLEAN_READS/"$READS_2" -k 21,31,51,71,81,101,127 -o $SCRATCH/2023_02_Sporothrix/04_assemblies/"$SAMPLE"
```

run Quast 5.0.2 on the assemblies for QC
20230605

```
#!/bin/bash

#SBATCH --account=def-tspribi
#SBATCH --time=1-0:0
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --job-name=quast
#SBATCH --output=quast.script.logs.out
#SBATCH --mail-user=w6p9c9j6t9c6a2i6@spribillelabworkspace.slack.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

SAMPLENAME="$1" #ca202112

module load StdEnv/2020
module load gcc/9.3.0
module load quast/5.0.2

quast --fungus -t 32 $SCRATCH/2023_02_Sporothrix/04_assemblies/"$SAMPLENAME"/contigs.fasta -o $SCRATCH/2023_02_Sporothrix/04_assemblies/"$SAMPLENAME"/quast
```
