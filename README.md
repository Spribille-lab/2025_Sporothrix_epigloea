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
* *Tremella yokohamensis* CBS 18117
Environmental sample:
* *Tremella yokohamensis* sporocarp TF4-1

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

* For CBS 139747 I combined bins 1, 2, 3, 21, 23, 26, and 34
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

### 5.1 ran funannotate pipeline with script `scripts/funannotate_annotate_script_v2023_1.sh`

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
cd "$project_directory"/10_annotations/DeepLoc2
conda activate deeploc-2.0
for sample in $(cat samples.txt)
do
deeploc2 -f ../protein_fasta/"$sample".proteins.fa -o "$sample" --model Accurate &>> "$sample".out
done
```

### 5.3 starfish v1.0.0

## 6. Genome comparisons
* compare utility of funannotate
* custom R scripts
    * ComplexHeatmap v2.13.2 package (Gu, 2022; Gu et al., 2016)
    * ape v5.6-2 (Paradis & Schliep, 2019)
    * DECIPHER v2.25.2 (Wright, 2016)
    * tidyverse v1.3.2 (Wickham et al., 2019)
* OrthoFinder 2.5.5 (Emms & Kelly, 2019)
* UpSetR v1.4.0 (Conway et al., 2017)
* IQ-TREE 2.0.7 (Minh et al., 2020)
    * ModelFinder (Kalyaanamoorthy et al., 2017)
    * 1000 ultrafast bootstrap replicates (Hoang et al., 2017)

### 6.1 comparisons
### 6.2 orthogroup inference
### 6.3 tree (show MSA)

## gene trees
