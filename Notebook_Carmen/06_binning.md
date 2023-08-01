# binned contigs using concoct 1.1.0
* 2023.06.11
* debary /data/ccallen/2023_02_Sporothrix/05_mapping/.

* mapped with:
* bowtie2 version 2.4.4 (graham)
* SAMtools version 1.16.1 (graham)


#ran concoct 1.1.0 in debary
```
cd /data/ccallen/2023_02_Sporothrix/06_binning
for sample in $(cat ../samples.txt)
do
bash ../scripts/06_run_concoct_on_single_metagenome.sh 2023_02_Sporothrix "$sample" &>> "$sample"_concoct_binning_script.out
echo just finished "$sample"
done
```

run_concoct_on_single_metagenome.sh

```
#!/bin/bash

project_directory="$1"
sample="$2"
#reads_1=/data/ccallen/"$project_directory"/03_CLEAN_READS/"$sample"_1.fastq.gz
#reads_2=/data/ccallen/"$project_directory"/03_CLEAN_READS/"$sample"_2.fastq.gz
contigs=/data/ccallen/"$project_directory"/04_assemblies/"$sample"/contigs.fasta
mapping_dir=/data/ccallen/"$project_directory"/05_mapping/"$sample"
binning_dir=/data/ccallen/"$project_directory"/06_binning/"$sample"

#bowtie2 version in narval is 2.4.4
#samtools version in narval is 1.16.1

eval "$(conda shell.bash hook)"
conda deactivate
eval "$(conda shell.bash hook)"
conda activate concoct-1.1.0

# I had already completed the mapping in compute canada so I copied those bam files to 05_mapping

# cut contigs into smaller parts
mkdir "$binning_dir"
cd "$binning_dir"
cut_up_fasta.py "$contigs" -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa

# Generate table with coverage depth information per sample and subcontig. This step assumes the directory
# ‘mapping’ contains sorted and indexed bam files where each sample has been mapped against the original contigs

concoct_coverage_table.py contigs_10K.bed "$mapping_dir"/"$sample".bam > coverage_table.tsv

# run concoct
concoct --composition_file contigs_10K.fa --coverage_file coverage_table.tsv -t 28 -b concoct_output/

# Merge subcontig clustering into original contig clustering
merge_cutup_clustering.py concoct_output/clustering_gt1000.csv > concoct_output/clustering_gt1000_merged.csv

# Extract bins as individual FASTA
mkdir concoct_output/fasta_bins
extract_fasta_bins.py "$contigs" concoct_output/clustering_gt1000_merged.csv --output_path concoct_output/fasta_bins
```