# Massive gene loss in the fungus *Sporothrix epigloea* accompanied a shift to life in a glucuronoxylomannan-based gel matrix

This repository contains scripts and intermediate results for the manuscript (Allen et al. 2025, Genome Biology and Evolution)

## Abstract
Fungi are well known for their ability to both produce and catabolize complex carbohydrates to acquire carbon, often in the most extreme of environments. Glucuronoxylomannan (GXM)-based gel matrices are widely produced by fungi in nature and though they are of key interest in medicine and pharmaceuticals, their biodegradation is poorly understood. Though some organisms, including other fungi, are adapted to life in and on GXM-like matrices in nature, they are almost entirely unstudied, and it is unknown if they are involved in matrix degradation. *Sporothrix epigloea* is an ascomycete fungus that completes its life cycle entirely in the short-lived secreted polysaccharide matrix of a white jelly fungus, *Tremella fuciformis*. To gain insight into how *S. epigloea* adapted to life in this unusual microhabitat, we compared the predicted protein composition of *S. epigloea* to that of 21 other *Sporothrix* species. We found that the genome of *S. epigloea* is smaller than that of any other sampled *Sporothrix*, with widespread functional gene loss, including those coding for serine proteases and biotin synthesis. In addition, many predicted CAZymes degrading both plant and fungal cell wall components were lost while a lytic polysaccharide monooxygenase (LPMO) with no previously established activity or substrate specificity, appears to have been gained. Phenotype assays suggest narrow use of mannans and other oligosaccharides as carbon sources. Taken together, the results suggest a streamlined machinery, including potential carbon sourcing from GXM building blocks, facilitates the hyperspecialized ecology of *S. epigloea* in the GXM-like milieu.
## Overview
```
project
├── README.md			# this doc; description of the repo and the project log
├── scripts			# all scripts generated for the analysis, with the exception of snakemake pipelines
├── 04_assemblies		
├── 08_GC_coverage		
├── 09_orthology
├── 10_annotations		
├── 11_LPMO_trees
├── 12_starfish
├── 14_CAFE 		
└── Notebook 			# temp log files for various parts of the analysis; the cleaned-up version of the same logs is in this file
└── results 			# results included in the manuscript
    ├── figures 		# figures
    └── tables 			# tables
```

## 1. Short-read sequencing, assembly, and binning

Isolates:
* *Sporothrix epigloea* CBS 119000
* *Sporothrix epigloea* CBS 573.63
* *Sporothrix bragantina* CBS 474.91
* *Sporothrix curviconia* CBS 959.73
* *Sporothrix eucalyptigena* CBS 139899
* *Sporothrix eucalyptigena* CBS 140593
* *Sporothrix thermara* CBS 139747
* *Annulohypoxylon annulatum* CBS 149473
* *Tremella fuciformis s.lat.* CBS 18117

Environmental sample:
* *Tremella fuciformis* sporocarp TF4-1

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
metawrap read_qc \
    -t 28 \
    -1 01_RAW_READS/"$sample"_1.fastq \
    -2 01_RAW_READS/"$sample"_2.fastq \
    -o 02_READ_QC/"$sample"
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

* assembly statistics generated with script `scripts/04_quast.sh`
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

## 2. Taxonomic assignments of bins
Software used:
* BUSCO 5.3.2 (Manni et al., 2021) - *Tremella* sporocarp
* BUSCO 5.4.7 (Manni et al., 2021) - CBS 139747

```
cd "$project_directory"/07_bin_QC
for path in $(ls bins/*.fa)
do
file=${path##*/}
bin=${file%.*}
echo $bin  
docker run -u $(id -u) -v $(pwd):/busco_wd ezlabgva/busco:v5.4.7_cv1 \
    busco \
    -i "$path" \
    -o "$bin" \
    -m genome \
    --cpu 8 \
    -f
done

# compile report
for filename in */*.txt
do
echo -n  $filename >> busco_report.txt;
sed -n -e '/^\tC:/p' $filename >> busco_report.txt
done
```

For CBS 139747 I combined bins 1, 2, 3, 21, 23, 26, and 34
    * C:95.6%[S:95.4%,D:0.2%],F:1.2%,M:3.2%,n:3817 (sordariomycetes_odb10)

For TF4-1 I extracted bin 63
    * C:98.0%[S:98.0%,D:0.0%],F:1.6%,M:0.4%,n:255 (fungi_odb10)

## 3. Generate GC and coverage plots
Software used:
* R libraries
    * tidyverse
    * Biostrings
    * DT
    * plotly

Used script `scripts/08_GC_cov.R` to generate GC plots.  Script was modified from Tagirdzhanova et al. (2021).

Plots can be viewed in `08_GC_coverage`


## 4. Gene prediction
Software used:
* funannotate v1.8.5
    * GeneMark-ES v4.68 for fungal genomes (Ter-Hovhannisyan et al., 2008)
    * AUGUSTUS v3.3.3 (Stanke et al., 2008)
    * SNAP v2006-07-28 (Korf, 2004)
    * GlimmerHMM v3.0.4 (Majoros et al., 2004) trained on BUSCO gene models
    * EVidenceModeler v1.1.1 (Haas et al., 2008)

ran with script `scripts/funannotate_predict_script_v2022_2.sh`

```
FILENAME="${1##*/}" #Tremella_mesenterica_S_ATCC28783_GCA_004117975.fna
DIRNAME=${FILENAME%.*} #Tremella_mesenterica_S_ATCC28783_GCA_004117975
SP=${FILENAME%%_S_*} #Tremella_mesenterica
SPECIES=${SP//_/ } #Tremella mesenterica
G=${SP%%_*} #Tremella
GS=${G:0:1} #T
S=${SP#*_} #mesenterica
ss=${S:0:2} #me
SS=$(tr '[a-z]' '[A-Z]' <<< $ss) #ME
SP_STRAIN=${FILENAME%%_G*} #Tremella_mesenterica_S_ATCC28783
STRAIN=${SP_STRAIN##*_} #ATCC28783
BASENAME="$GS""$SS""$STRAIN" #TMEATCC28783

eval "$(conda shell.bash hook)"
conda deactivate
eval "$(conda shell.bash hook)"
conda activate fun-env #version 1.8.5
export FUNANNOTATE_DB=/data/databases/funannotate_db
export PATH=$PATH:/home/ccallen/.local/bin/signalp-4.1 #Change for your path to signalp
mkdir "$DIRNAME"
funannotate clean -i ../00_RawData/"$FILENAME" -o  "$DIRNAME"/"$DIRNAME"_cleaned.fas --exhaustive &&
funannotate sort -i  "$DIRNAME"/"$DIRNAME"_cleaned.fas -o  "$DIRNAME"/"$DIRNAME"_sorted.fas &&
funannotate mask -i  "$DIRNAME"/"$DIRNAME"_sorted.fas -o  "$DIRNAME"/"$DIRNAME"_masked.fas --cpus 12
mv funannotate-mask.log ./"$DIRNAME"/
mkdir ~/gmes_tmp/"$DIRNAME"_gmes #Change if needed
cp "$DIRNAME"/"$DIRNAME"_masked.fas ~/gmes_tmp/"$DIRNAME"_gmes/
cd ~/gmes_tmp/"$DIRNAME"_gmes
~/.local/bin/gmes_linux_64/gmes_petap.pl --ES --fungus --cores 12 --sequence "$DIRNAME"_masked.fas #Change gmes_petap.pl path if needed
cd /data/ccallen/2022_04_fungal_annotation/01_funannotate #Change for your project path
mv ~/gmes_tmp/"$DIRNAME"_gmes/ ./"$DIRNAME"
funannotate predict -i "$DIRNAME"/"$DIRNAME"_masked.fas -s "$SPECIES" -o "$DIRNAME"/"$BASENAME"_preds --optimize_augustus --cpus 12 --genemark_gtf "$DIRNAME"/"$DIRNAME"_gmes/genemark.gtf --name "$BASENAME" #if evidence modeler fails you can try adding --no-evm-partitions
```

* Estimated genome completion using EukCC2 run on proteins predicted by funannotate

```
cd 10_annotations
conda activate eukcc-2.1.0
for file in 10_annotations/protein_fasta/*.fa
do
fasta=${file##*/}
genome=${fasta%.proteins.fa}
mkdir eukcc2/"$genome"
eukcc single --out eukcc2/"$genome" --threads 8 protein_fasta/"$fasta" --db /data/databases/eukcc2_db/eukcc2_db_ver_1.1
done
```

## 5. Genome annotation
Software used:
* funannotate v1.8.15
    * InterProScan5.63-95.0 (Blum et al., 2021)
    * eggNOG-mapper v2.1.11 (Cantalapiedra et al., 2021)
    * Pfam 35.0 R package version 2.6-2 (El-Gebali et al., 2019)
    * UniProt DB v2023_03 (The UniProt Consortium, 2019)
    * dbCAN 12.0 (Zheng et al., 2023)
    * MEROPS v12.0 (Rawlings et al., 2018)
    * BUSCO v2.0.0 dikarya models (Dykaria_odb9; Simão et al., 2015)
    * Gene2Product database v1.91 (https://github.com/nextgenusfs/gene2product)
    * antiSMASH v7 (fungal version; Blin et al., 2023)
    * SignalP v4.1 (Petersen et al., 2011)
* DeepLoc 2.0 (Thumuluri et al., 2022)
* starfish v1.0.0 (Gluck-Thaler & Vogan, 2023)
    * metaeuk (Karin et al., 2020)
* kofamscan v1.3.0 (Aramaki et al., 2019)
* ggkegg 1.1.18 (Sato et al., 2023)

### 5.1 ran funannotate pipeline
script `scripts/funannotate_annotate_script_v2023_1.sh`

```
FILENAME="${1##*/}" #Tremella_mesenterica_S_ATCC28783_GCA_004117975.fna
DIRNAME=${FILENAME%.*} #Tremella_mesenterica_S_ATCC28783_GCA_004117975
SP=${FILENAME%%_S_*} #Tremella_mesenterica
SPECIES=${SP//_/ } #Tremella mesenterica
G=${SP%%_*} #Tremella
GS=${G:0:1} #T
S=${SP#*_} #mesenterica
ss=${S:0:2} #me
SS=$(tr '[a-z]' '[A-Z]' <<< $ss) #ME
SP_STRAIN=${FILENAME%%_G*} #Tremella_mesenterica_S_ATCC28783
STRAIN=${SP_STRAIN##*_} #ATCC28783
BASENAME="$GS""$SS""$STRAIN" #TMEATCC28783

# set the assembly directory (input) and the funannotate directory (output)
ASSEMBLY_DIR="/data/ccallen/assemblies"
FUNANNOTATE_DIR="/data/ccallen/2023_02_Sporothrix/10_annotations/funannotate"

# export paths so that funannotate can find these things
export FUNANNOTATE_DB=/data/databases/funannotate_db_2023
export PATH=$PATH:/home/ccallen/.local/bin/signalp-4.1 #Change for your path to signalp
export PATH=$PATH:/data/ccallen/bin/phobius #Change to your path

# Set conda environment to funannotate
eval "$(conda shell.bash hook)"
conda deactivate
eval "$(conda shell.bash hook)"
conda activate funannotate-1.8.15-env #version 1.8.15

# change into your funannotate directory and make a directory named the same as the assembly name
cd "$FUNANNOTATE_DIR"
#mkdir "$DIRNAME"

# run funannotate initial steps (clean, sort, mask) to setup for predictions
#funannotate clean -i "$ASSEMBLY_DIR"/"$FILENAME" -o "$DIRNAME"/"$DIRNAME"_cleaned.fas --exhaustive &&
#funannotate sort -i  "$DIRNAME"/"$DIRNAME"_cleaned.fas -o "$DIRNAME"/"$DIRNAME"_sorted.fas &&
#funannotate mask -i  "$DIRNAME"/"$DIRNAME"_sorted.fas -o "$DIRNAME"/"$DIRNAME"_masked.fas --cpus 12
#mv funannotate-mask.log ./"$DIRNAME"/


# run Genemark in the home directory and move the output to the funannotate directory
#mkdir ~/gmes_tmp/"$DIRNAME"_gmes #Change if needed
#cp "$DIRNAME"/"$DIRNAME"_masked.fas ~/gmes_tmp/"$DIRNAME"_gmes/
#cd ~/gmes_tmp/"$DIRNAME"_gmes
#~/.local/bin/gmes_linux_64/gmes_petap.pl --ES --fungus --cores 12 --sequence "$DIRNAME"_masked.fas #Change gmes_petap.pl path if needed
#cd "$FUNANNOTATE_DIR"
#mv ~/gmes_tmp/"$DIRNAME"_gmes/ ./"$DIRNAME"

# predict gene models
#funannotate predict -i "$DIRNAME"/"$DIRNAME"_masked.fas -s "$SPECIES" -o "$DIRNAME"/"$BASENAME"_preds --optimize_augustus --cpus 12 --genemark_gtf "$DIRNAME"/"$DIRNAME"_gmes/genemark.gtf --name "$BASENAME" #if evidence modeler fails you can try adding --no-evm-partitions

# run iprscan in funannotate
funannotate iprscan -i "$DIRNAME"/"$BASENAME"_preds -m local --iprscan_path /data/databases/interproscan-5.63-95.0/interproscan.sh -c 12 &&

# run antismash
eval "$(conda shell.bash hook)"
conda deactivate
eval "$(conda shell.bash hook)"
conda activate antismash-latest #6.1.1 atm
antismash --taxon fungi --output-dir "$DIRNAME"/"$BASENAME"_preds/annotate_misc/antismash --genefinding-tool none "$DIRNAME"/"$BASENAME"_preds/predict_results/"$SP".gbk &&

# run eggnog mapper with version 5 database
eval "$(conda shell.bash hook)"
conda deactivate
eval "$(conda shell.bash hook)"
conda activate eggnog-mapper-2.1.11-env #Change to your own emapper envionment
emapper.py -i "$DIRNAME"/"$BASENAME"_preds/predict_results/"$SP".proteins.fa --output "$DIRNAME"/"$BASENAME"_preds/annotate_misc/eggnog_results -d euk --data_dir /data/databases/eggnogdb5 --cpu 12 --override -m diamond &&

# run annotate
eval "$(conda shell.bash hook)"
conda deactivate
eval "$(conda shell.bash hook)"
conda activate funannotate-1.8.15-env
#phobius.pl -short "$DIRNAME"/"$BASENAME"_preds/predict_results/"$SP".proteins.fa > "$DIRNAME"/"$BASENAME"_preds/annotate_misc/phobius.results.txt
funannotate annotate -i "$DIRNAME"/"$BASENAME"_preds --sbt template.sbt.txt --eggnog "$DIRNAME"/"$BASENAME"_preds/annotate_misc/eggnog_results.emapper.annotations --antismash "$DIRNAME"/"$BASENAME"_preds/annotate_misc/antismash/"$SP".gbk --busco_db dikarya --strain "$STRAIN" --cpus 12 #Change location of template file if needed.  You may want to choose a bifferent --busco_db. You may not need the --strain flag.

#clean-up
mv iprscan_* "$DIRNAME"/.
```
* ran GBK files through antiSMASH 7 fungi version on the webserver

### 5.2 Deeploc

```
cd 10_annotations/DeepLoc2
conda activate deeploc-2.0
for sample in $(cat samples.txt)
do
deeploc2 -f ../protein_fasta/"$sample".proteins.fa -o "$sample" --model Accurate &>> "$sample".out
done
```

### 5.3 starfish

```
cd 12_starfish
```
* create ome2*.txt files detailing the absolute path to each genome's gff3 and assembly
```
realpath assembly/*.final.fasta | perl -pe 's/^(.+?([^\/]+?).contigs.final.fasta)$/\2\t\1/' > ome2assembly.txt
realpath gff3/*.final.gff3 | perl -pe 's/^(.+?([^\/]+?).final.gff3)$/\2\t\1/' > ome2gff.txt
```

* concatenate all assembly files and make a blastn database:
```
mkdir blastdb
cut -f2 ome2assembly.txt | xargs cat > blastdb/spor29.assemblies.fna
makeblastdb \
-in blastdb/spor29.assemblies.fna \
-out blastdb/spor29.assemblies \
-parse_seqids \
-dbtype nucl
```

* *de novo* annotate tyrs with the provided YR HMM and amino acid queries
```
mkdir geneFinder
starfish annotate \
-T 28 \
-x spor29_tyr \
-a ome2assembly.txt \
-g ome2gff.txt \
-p /data/ccallen/bin/starfish/database/YRsuperfams.p1-512.hmm \
-P /data/ccallen/bin/starfish/database/YRsuperfamRefs.faa \
-i tyr \
-o geneFinder/
```

* results in `12_starfish`

### 5.4 KEGG ortholog assignment

* get the list of KEGG modules
```
cd /Users/carmenallen/Documents/2023_02_Sporothrix/10_annotations/kegg
wget https://rest.kegg.jp/list/module
```

* get module files from the KEGG API
```
cd mf
while IFS=$'\t' read -r v1 v2; do wget https://rest.kegg.jp/get/"$v1"; sleep 1; done < ../module
```

* run kofamscan
```
/data/ccallen/bin/kofam_scan/exec_annotation -o kegg_relaxed/{sample}/{sample}.ko.txt -f detail-tsv --threshold-scale=0.75 --tmp-dir tmp_d_{wildcards.sample} protein_fasta/{sample}.proteins.fa
/data/ccallen/bin/kofam_scan/exec_annotation -o kegg_relaxed/{sample}/{sample}.kegg.mapper.txt -f mapper --threshold-scale=0.75 --tmp-dir tmp_m_{wildcards.sample} protein_fasta/{sample}.proteins.fa
```
* counted the number of KOs per genome using `scripts/10_KEGG_count_identifiers.R`
* calculate kegg module completeness using ggkegg, create heatmap and completeness figure `scripts/10_ggkegg_with_partners_v2.R`.  Script was modified from Tagirdzhanova et al. (2024).
* results in `10_annotations/kegg_relaxed`

### 5.5 merged annotations from funannotate, DeepLoc, Orthofinder, and KEGG

* merged annotations using `scripts/10_merge_annotations.R`

## 6. Genome comparisons
Software used:
* compare utility of funannotate
* OrthoFinder 2.5.5 (Emms & Kelly, 2019)
* KinFin 1.1.1 (Laetsch & Blaxter, 2017)
* CAFE5 1.1 (Mendes et al., 2020)
* CafePlotter (https://github.com/moshi4/CafePlotter)
* custom R scripts
    * ComplexHeatmap v2.13.2 package (Gu, 2022; Gu et al., 2016)
    * ape v5.6-2 (Paradis & Schliep, 2019)
    * DECIPHER v2.25.2 (Wright, 2016)
    * tidyverse v1.3.2 (Wickham et al., 2019)
    * UpSetR v1.4.0 (Conway et al., 2017)
    * phylotools v0.2.2 (Zhang et al.)
    * ggtree 3.10.1 (Xu et al., 2022)
    * phytools 2.3-0 (Revell, 2024)
    * treeio 1.26.0 (Wang et al., 2020)
    * tidytree 0.4.6 (Yu, 2022)
    * ggVennDiagram (Gao & Dusa, 2024)
    * cowplot v1.1.3 (Wilke et al., 2024)
* IQ-TREE 2.0.7 (Minh et al., 2020)
    * ModelFinder (Kalyaanamoorthy et al., 2017)
    * ultrafast bootstraps (Hoang et al., 2017)

### 6.1 Species tree construction of *Sporothrix*, *Ophiostoma*, and *Leptographium* genomes

* orthogroup inference `scripts/09_orthofinder_docker_for_tree.sh`
* results found in `09_orthology/OrthoFinder/Results_Aug03/`
```
orthofinder \
-t 64 \
-a 64 \
-M msa \
-A mafft \
-o /input/orthofinder_output_Aug3 \
-f /input/orthofinder_input
```

* run IQ-TREE on MSA from orthofinder output using `scripts/09_iqtreesh.sh`
* MSA input can be found in `09_orthology/OrthoFinder/Results_Aug03/MultipleSequenceAlignments/`
* results can be found in `09_orthology/OrthoFinder/Results_Aug03/Species_Tree/iqtree2`

```
iqtree2 -s Results_Aug03/MultipleSequenceAlignments/SpeciesTreeAlignment.fa \
-o Leptographium_lundbergii_CBS138716 \
-nt 29 \
-m MFP \
-B 1000 \
-pre sporothrix_placement
```

### 6.2 Orthogroup inference of *Sporothrix* *sensu stricto* genomes
* using `scripts/09_orthofinder_core_docker.sh`
* results found in `09_orthology/OrthoFinder/Results_Aug07/`

```
orthofinder \
-t 64 \
-a 64 \
-M msa \
-A mafft \
-o /input/orthofinder_output_core_Aug7 \
-f /input/orthofinder_input_core
```

* used KinFin to produce representative annotation per orthogroup.

```
cd "$project_dir"/09_orthology/orthofinder_input_core
cat *.fa > all_Sporothrix_proteins.fa

#run InterProScan
conda activate funannotate-1.8.15-env #to get java to work
cd "$project_dir"/09_orthology/KinFin
/data/databases/interproscan-5.63-95.0/interproscan.sh -i "$project_dir"/09_orthology/orthofinder_input_core/all_Sporothrix_proteins.fa -d out/ -dp -t p --goterms -appl Pfam -f TSV
conda deactivate

#covert to table readable by KinFin
conda activate kinfin
kinfin/scripts/iprs2table.py -i out/all_Sporothrix_proteins.fa.tsv

#run KinFin
kinfin -g Orthogroups.txt -c config.txt -s SequenceIDs.txt -p SpeciesIDs.txt -f functional_annotation.txt -a protein_fasta/
kinfin/scripts/functional_annotation_of_clusters.py all -f kinfin_results/cluster_domain_annotation.IPR.txt -c kinfin_results/cluster_counts_by_taxon.txt --domain_protein_cov 0.3 --domain_taxon_cov 0.02 -o IPR
kinfin/scripts/functional_annotation_of_clusters.py all -f kinfin_results/cluster_domain_annotation.GO.txt -c kinfin_results/cluster_counts_by_taxon.txt --domain_protein_cov 0.3 --domain_taxon_cov 0.02 -o GO
kinfin/scripts/functional_annotation_of_clusters.py all -f kinfin_results/cluster_domain_annotation.Pfam.txt -c kinfin_results/cluster_counts_by_taxon.txt --domain_protein_cov 0.3 --domain_taxon_cov 0.02 -o Pfam
```

* visualized the orthogroup overlaps using `scripts/09_Sporothrix_orthofinder_KinFin_results.R`
* extracted orthogroups unique to *S. epigloea* and undetected in *S. epigloea* using `scripts/09_Sporothrix_orthofinder_KinFin_results.R`
* manually categorized the InterPro annotations and added a KEGG annotation if at least 2 (for those unique to *S. epigloea*) or 5 proteins (for those not found in *S. epigloea*) in the orthogroup contained that annotation  `09_orthology/KinFin/OrthoFinder/combined_cluster_functional_annotation.gained.p30.x2.xls` and `09_orthology/KinFin/OrthoFinder/combined_cluster_functional_annotation.lost.p30.x2.xls`
* visualized the orthogroup annotations using `scripts/09_Sporothrix_orthofinder_KinFin_results.R`

### 6.3 Comparisons using funannotate

* used compare utility of funannotate to compare annotations
```
funannotate compare -i \
<list gbk files> \
--num_orthos 1 --cpus 32 -o funannotate_compare_notree
```

### 6.4 visualize annotations

* visualized the genome comparisons using `scripts/10_comparative_genomic_analysis.R`

### 6.5 Gene family expansion and contraction analysis

* Coerced phylogenomic tree into ultrametric tree using `scripts/14_CAFE.R`
* trimmed phylogenomic tree for orthogroups using `scripts/14_CAFE.R`
* prepared CAZyme data for CAFE using `scripts/10_process_dbcanV12.R`
* prepared protease, BGCs, and orthogroup data for CAFE using `scripts/14_CAFE.R`


* test for the number of discrete rate categories (*K*) that produces the greatest likelihood
```
category=cazymefamily

tree="$project_dir"/09_orthology/Results_Aug03/Species_Tree/iqtree2/sporothrix_placement.rerooted.ultrametric.contree
input="$project_dir"/10_annotations/dbcan4/summarized_outputs/cazy_summarized_CAFE.txt

# run CAFE5
for n in {1..10}
do
output=cazymefamily_singlelambda_k"$n"
mkdir "$output"
/data/ccallen/bin/CAFE5/bin/cafe5 -t "$tree" -i "$input" -p -k "$n" -o "$output" &>> "$output"/cafe.log
done
```

* compile -lnL for each into one document
```
for n in {1..10}
do
echo "$category"_singlelambda_k"$n" >> "$category"_discrete_rate_categories.txt
grep "Model" "$category"_singlelambda_k"$n"/Base_results.txt >> "$category"_discrete_rate_categories.txt
grep "Model" "$category"_singlelambda_k"$n"/Gamma_results.txt >> "$category"_discrete_rate_categories.txt
done
```

* count significant families at the p=0.05 threshold
```
for n in {1..10}
do
echo "$category"_singlelambda_k"$n" >> "$category"_significant_families_count.txt
grep -c "y" "$category"_singlelambda_k"$n"/*_family_results.txt >> "$category"_significant_families_count.txt
done
```

* results are found in `14_CAFE/`

Due to failure of convergence with the gamma models the Base model was used (*K* = 1)

CAZyme families, Model Base Final Likelihood (-lnL): 4159.66, Lambda: 0.92955594638637
Protease families, Model Base Final Likelihood (-lnL): 2676.62, Lambda: 0.60118488521017
BGCs, Model Base Final Likelihood (-lnL): 280.253, Lambda: 0.99761242999222
Orthogroups, Model Base Final Likelihood (-lnL): 115429, Lambda: 0.57113146703192

* visualized using CafePlotter

```
cafeplotter -i $cafe_folder -o $output --format png --fig_width 20 --fig_height 0.4 --count_label_size 10
```
* also visualized the results with a custom R script found at `scripts/14_CAFE.R`
* the significant expansion and contraction events at the *S. epigloea* tips and two transition nodes were parsed in `scripts/14_CAFE.R` and presented in Supplementary Tables S5 and S6.
* annotations for significant contraction events were extracted from annotation tables using `scripts/09_Sporothrix_orthofinder_KinFin_results.R`

## 7. LPMO gene trees
Software used:
* dbCAN 4.0.0 (Zheng et al., 2023)
* MAFFT v7.450 (Katoh & Standley, 2013)
* trimal v1.4.1 (Capella-Gutierrez et al., 2009)
* IQ-TREE 2.0 (Minh et al., 2020)

* used dbCAN 4.0.0 (dbCAN HMMdb release 12.0) to annotate CAZymes from 10 *de novo* genome assemblies and 58 additional downloaded assemblies 

Agaricus_bisporus_S_H97_GCA_000300575
Ampelomyces_quisqualis_S_BRIP72107_GCA_018398575
Annulohypoxylon_annulatum_S_SEMC78_GCA_944038255
Annulohypoxylon_stygium_S_MG137_GCA_003049155
Annulohypoxylon_truncatum_S_CBS140778_GCA_902805465
Arthrobotrys_oligospora_S_ATCC24927_GCA_000225545
Cladobotryum_protrusum_S_CCMJ2080_GCA_004303015
Clonostachys_rosea_S_671_GCA_000963775
Cryptococcus_gattii_S_WM276_GCA_000185945
Cryptococcus_neoformans_S_V31_GCA_002215745
Escovopsis_weberi_S_strain_GCA_001278495
Holtermannia_corniformis_S_JCM1743_GCA_001599935
Hypomyces_perniciosus_S_HP10_GCA_008477525
Leptographium_lundbergii_S_CBS138716_GCA_001455505
Ophiostoma_fasciatum_S_VPRI43845_GCA_019925495
Ophiostoma_ips_S_VPRI43529_GCA_019925475
Ophiostoma_novoulmi_S_H327_GCA_000317715
Parasitella_parasitica_S_CBS41266_GCA_000938895
Pseudozyma_flocculosa_S_PF1_GCA_000417875
Pycnoporus_coccineus_S_BRFM310_GCA_002092935
Sporothrix_bragantina_S_CBS47491_GCA_XXXXXXXXX
Sporothrix_brasiliensis_S_5110_GCA_000820605
Sporothrix_brunneoviolacea_S_CBS124561_GCA_021396205
Sporothrix_curviconia_S_CBS95973_GCA_XXXXXXXXX
Sporothrix_dimorphospora_S_CBS55374_GCA_021397985
Sporothrix_epigloea_S_CBS119000_GCA_943908295
Sporothrix_epigloea_S_CBS57363_GCA_943900835
Sporothrix_epigloea_S_TF4163_GCA_944036445
Sporothrix_eucalyptigena_S_CBS139899_GCA_XXXXXXXXX
Sporothrix_eucalyptigena_S_CBS140593_GCA_XXXXXXXXX
Sporothrix_euskadiensis_S_VPRI43754_GCA_019925375
Sporothrix_globosa_S_CBS120340_GCA_001630435
Sporothrix_humicola_S_CBS118129_GCA_021396245
Sporothrix_inflata_S_CBS23968_GCA_021396225
Sporothrix_insectorum_S_RCEF264_GCA_001636815
Sporothrix_luriei_S_CBS93772_GCA_021398005
Sporothrix_mexicana_S_CBS120341_GCA_021396375
Sporothrix_nigrograna_S_VPRI43755_GCA_019925305
Sporothrix_pallida_S_CBS13156_GCA_021396235
Sporothrix_phasma_S_CBS119721_GCA_011037845
Sporothrix_protearum_S_CBS116654_GCA_016097115
Sporothrix_pseudoabietina_S_VPRI43531_GCA_019925295
Sporothrix_schenckii_S_1099_GCA_000961545
Sporothrix_thermara_S_CBS139747_GCA_XXXXXXXXX
Sporothrix_variecibatus_S_CBS121961_GCA_021396255
Tremella_fuciformis_S_NBRC9317_GCA_024343405
Tremella_fuciformis_S_SEMC49_GCA_944472115
Tremella_fuciformis_S_tr26_GCA_000987905
Tremella_mesenterica_S_ATCC28783_GCA_004117975
Tremella_mesenterica_S_DSM1558_GCA_000271645
Tremella_yokohamensis_S_NBRC100148_GCA_024343385
Trichoderma_atroviride_S_IMI206040_GCA_000171015
Trichoderma_harzianum_S_CBS22695_GCA_003025095
Trichoderma_parareesei_S_CBS125925_GCA_001050175
Trichoderma_reesei_S_QM6a_GCA_000167675
Trichoderma_virens_S_Gv298_GCA_000170995
Kockovaella_imperatae_S_NRRLY17943_GCA_002102565
Melampsora_medusae_S_Mmd05TRE539_GCA_002157035
Ophiocordyceps_australis_S_Map64_GCA_002591415
Rhizopus_arrhizus_S_971192_GCA_000697195
Saccharomyces_cerevisiae_S_S288C_GCA_000146045
Schizosaccharomyces_japonicus_S_yFS275_GCA_000149845
Sclerotinia_sclerotiorum_S_1980_GCA_000146945
Talaromyces_cellulolyticus_S_NA01_GCA_009805475
Tuber_melanosporum_S_Mel28_GCA_000151645
Ustilago_maydis_S_521_GCA_000328475
Xylaria_hypoxylon_S_CBS122620_GCA_902806585
Xylographa_opegraphella_S_Spribille41601_GCA_022813865

```
run_dbcan {input} protein \
--out_dir dbcan4/{wildcards.sample}_dbcan \
--db_dir /data/ccallen/bin/run_dbcan_4.0.0/db/ \
--use_signalP True
```
* processed the outputs and generated a table of lpmo and AA14 coordinates using `scripts/10_process_dbcan.R`

* extracted LPMO sequences from predicted proteomes

```
cd 10_annotations/dbcan4
while IFS=$'\t' read -r v1 v2 v3 v4; do samtools faidx ../protein_fasta_dbcan/"$v1".proteins.fa $v2:$v3-$v4 >> summarized_outputs/allgenomes.cazy_aa14.fa; done < summarized_outputs/AA14_coordinates.txt
while IFS=$'\t' read -r v1 v2 v3 v4; do samtools faidx ../protein_fasta_dbcan/"$v1".proteins.fa $v2:$v3-$v4 >> summarized_outputs/allgenomes.cazy_lpmo.fa; done < summarized_outputs/lpmo_coordinates.txt
```

### 7.1 LPMO tree

* aligned LPMO sequences with MAFFT v7.450 (within Geneious Prime® 2023.2.1) and trimmed with trimal v1.4.1

```
trimal -in LPMO_alignment.phy -out LPMO_alignment_trimmed.phy -gappyout
```

* generated a maximum likelihood tree of LPMOs with IQ-TREE 2.0-rc1 using ModelFinder and 1000 ultrafast bootstrap replicates

```
iqtree \
-s 11_LPMO_trees/LPMO/LPMO_alignment_trimmed.phy \
-m MFP \
-B 1000 \
-pre LPMO
```

* alignment and results in `11_LPMO_trees/LPMO`

### 7.2 AA14 tree

* used blastp to search the nr_clustered database
* accessed 2024.03.07
* query: original AA14 characterized AUM86167.1
* evalue cutoff of 1e-25
* 1308 blastp hits

* combined with 59 AA14 protein sequences from dbCAN 4.0.0 genome annotations
```
cat 11_LPMO_trees/AA14/Fungi_March2024/blastp/blastp_results.fa 10_annotations/dbcan4/summarized_outputs/allgenomes.cazy_aa14.renamed.fa > 11_LPMO_trees/AA14/Fungi_March2024/Fungi_AA14.fa
```

* aligned AA14 sequences with MAFFT v7.450 (within Geneious Prime® 2023.2.1) and trimmed with trimal v1.4.1

```
trimal -in Fungi_AA14_alignment.phy -out Fungi_AA14_alignment_trimmed.phy -gappyout
```

* generated a maximum likelihood tree of AA14s with IQ-TREE 2.0-rc1 using ModelFinder and 1000 ultrafast bootstrap replicates

```
cd 11_LPMO_trees/AA14/Fungi_March2024
iqtree \
-s Fungi_AA14_alignment_trimmed.phy \
-nt auto \
-m MFP \
-B 1000 \
-o PRP77934.1_Planoprotostelium_fungivorum,PRP85073.1_Planoprotostelium_fungivorum,PRP85061.1_Planoprotostelium_fungivorum,PRP85059.1_Planoprotostelium_fungivorum \
-pre AA14_Fungi_outgroup1
```

* alignment and results in `11_LPMO_trees/AA14/Fungi_March2024`
* used `scripts/11_taxize.R` to pull taxonomy for each tip and generate helper files for iTOL (Letunic & Bork 2019).
* text files used to annotate tree in `11_LPMO_trees/AA14/Fungi_March2024/itol20241022`

## 8. Phenotype arrays

* custom R scripts
    * tidyverse v1.3.2 (Wickham et al., 2019)
    * cowplot v1.1.3 (Wilke et al., 2024)

* visualized the orthogroup overlaps using `scripts/phenotype_arrays_Sporothrix.R`
