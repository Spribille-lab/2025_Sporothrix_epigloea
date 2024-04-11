# Gene loss in *Sporothrix epigloea* accompanied a shift to life in the *Tremella* fungal heteropolymer matrix

This repository contains scripts and intermediate results for the manuscript (Allen et al. 2024, <journal>)

## Abstract
The unique microhabitat of *Sporothrix epigloea*—the extracellular matrix of *Tremella* jelly fungi—is rich with acidic glucuronoxylomannan (GXM) polysaccharides.  GXM polysaccharides are medically relevant both as a pharmaceutical extracted from *Tremella* fruiting bodies and as a critical virulence factor involved in cryptococcosis caused by the tremelloid fungus *Cryptococcus*.  Using a comparative genomics approach, we surveyed the genome assemblies of *Sporothrix* species through the evolutionary lifestyle transition from soil or bark beetle associate to life on a GXM-rich fungal fruiting body.  In addition, we inferred protein orthogroups across core *Sporothrix* genomes to identify gene families that were either lost or derived in the divergence of *S. epigloea*.  We found that *S. epigloea* genomes were smaller than any other *Sporothrix* genome and reduced in carbohydrate active enzymes (CAZymes), proteases, biosynthetic gene clusters, and sugar transporters.  The suite of CAZyme families degrading both plant and fungal cell wall components were reduced in *S. epigloea*.  At the same time, a lytic polysaccharide monooxygenase (LPMO) with no previously established activity or substrate specificity, appears to have been uniquely acquired in *S. epigloea*.  This result calls for further investigation in the substrate activity of these enzymes and if they play a role in degrading GXM in *Tremella* fruiting bodies.


## Overview
```
project
├── README.md			# this doc; description of the repo and the project log
├── scripts				# all scripts generated for the analysis, with the exception of snakemake pipelines
├── 04_assemblies		
├── 08_GC_coverage		
├── 09_orthology
├── 10_annotations		
├── 11_LPMO_trees
├── 12_starfish			
└── Notebook 			# temp log files for various parts of the analysis; the cleaned-up version of the same logs is in this file
└── results 			# results included in the manuscript
    ├── figures 		# figures
    └── tables 			# tables
```

## 1. Assembly and binning
Software used:
* metaWRAP pipeline v.1.3.2 (Uritskiy et al. 2018)
* metaSPAdes v3.15.x (Nurk et al. 2017)
* QUAST v5.0.2 (Mikheenko et al., 2018)
* Bowtie v2.4.4 (Langmead & Salzberg, 2012)
* SAMtools v1.16.1 (Danecek et al. 2021)
* CONCOCT v1.1.0 (Alneberg et al. 2014)

### 1.1. Adapter and quality trimming
```
for sample in $(cat samples.txt)
do
metawrap read_qc -t 28 -1 01_RAW_READS/"$sample"_1.fastq -2 01_RAW_READS/"$sample"_2.fastq -o 02_READ_QC/"$sample"
done
```
### 1.2. Assembly

* used the script `scripts/04_spades_cc.sh`
```
project_dir=project/dir/path
SAMPLE="$1"
READS_1="$SAMPLE"_1.fastq.gz
READS_2="$SAMPLE"_2.fastq.gz

module load StdEnv/2020
module load spades/3.15.4

spades.py -t ${SLURM_CPUS_PER_TASK} \
    -m 2000 \
    --meta \
    --pe1-1 "$project_dir"/03_CLEAN_READS/"$READS_1" \
    --pe1-2 "$project_dir"/03_CLEAN_READS/"$READS_2" \
    -k 21,31,51,71,81,101,127 \
    -o "$project_dir"/04_assemblies/"$SAMPLE"
```

* assembly statistics with script `scripts/04_quast.sh`

```
project_dir=project/dir/path
SAMPLENAME="$1"

module load StdEnv/2020
module load gcc/9.3.0
module load quast/5.0.2

quast \
    --fungus \
    -t 32 \
    "$project_dir"/04_assemblies/"$SAMPLENAME"/contigs.fasta \
    -o "$project_dir"/04_assemblies/"$SAMPLENAME"/quast
```

### 1.3.  Binning

* mapped reads to contigs using script `scripts/05_run_mapping_on_single_metagenome`
```
module load StdEnv/2020
module load gcc/9.3.0
module load bowtie2/2.4.4
module load samtools/1.16.1

project_directory="$1"
sample="$2"
reads_1=$SCRATCH/"$project_directory"/03_CLEAN_READS/"$sample"_1.fastq.gz
reads_2=$SCRATCH/"$project_directory"/03_CLEAN_READS/"$sample"_2.fastq.gz
contigs=$SCRATCH/"$project_directory"/04_assemblies/"$sample"/contigs.fasta
mapping_dir=$SCRATCH/"$project_directory"/05_mapping/"$sample"

mkdir "$mapping_dir"
cd "$mapping_dir"
bowtie2-build "$contigs" "$sample"_assembly
bowtie2 -x "$sample"_assembly -1 "$reads_1" -2 "$reads_2" --no-unal -p 32 -S "$sample".sam
samtools view -b -o "$sample"-raw.bam "$sample".sam
samtools sort -o "$sample".bam "$sample"-raw.bam
samtools index "$sample".bam
rm "$sample".sam
rm "$sample"-raw.bam
```

* binned contigs using script `scripts/06_binning`
```
project_directory="$1"
sample="$2"
contigs=/data/ccallen/"$project_directory"/04_assemblies/"$sample"/contigs.fasta
mapping_dir=/data/ccallen/"$project_directory"/05_mapping/"$sample"
binning_dir=/data/ccallen/"$project_directory"/06_binning/"$sample"

mkdir "$binning_dir"
cd "$binning_dir"
cut_up_fasta.py "$contigs" -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa
concoct_coverage_table.py contigs_10K.bed "$mapping_dir"/"$sample".bam > coverage_table.tsv
concoct --composition_file contigs_10K.fa --coverage_file coverage_table.tsv -t 28 -b concoct_output/
merge_cutup_clustering.py concoct_output/clustering_gt1000.csv > concoct_output/clustering_gt1000_merged.csv
mkdir concoct_output/fasta_bins
extract_fasta_bins.py "$contigs" concoct_output/clustering_gt1000_merged.csv --output_path concoct_output/fasta_bins
```