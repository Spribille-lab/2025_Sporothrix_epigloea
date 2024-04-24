

# annotate eukaryote genomes based on KEGG Orthology and hidden Markov models
* using KofamScan version 1.3.0
* 2024.04.23

to download KofamScan version 1.3.0
```
cd /data/databases/KOfam_2024
wget ftp://ftp.genome.jp/pub/db/kofam/*
gzip -d ko_list.gz
tar -xvzf profiles.tar.gz
cd /data/ccallen/bin/
git clone "https://github.com/takaram/kofam_scan"
-> put profile: /data/databases/KOfam_2024/profiles/eukaryote in the config.yml
```

```
cd /data/ccallen/2023_02_Sporothrix/10_annotations
cat genomes.txt

genome
Ophiostoma_fasciatum_S_VPRI43845_GCA_019925495
Sporothrix_brasiliensis_S_5110_GCA_000820605
Sporothrix_bragantina_S_CBS47491_GCA_XXXXXXXXX
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
Leptographium_lundbergii_S_CBS138716_GCA_001455505
Ophiostoma_ips_S_VPRI43529_GCA_019925475
Ophiostoma_novoulmi_S_H327_GCA_000317715
```

# snakefile 
* /data/ccallen/2023_02_Sporothrix/10_annotations
```
import glob
import os
import pandas as pd

genome_ann = pd.read_table('genomes.txt').set_index("genome", drop=False)
GENOMES_ANN = list(genome_ann['genome'])

genome_dbcan = pd.read_table('genomes_dbcan.txt').set_index("genome", drop=False)
GENOMES_DBCAN = list(genome_dbcan['genome'])

rule all:
        input:
                expand("kegg/{sample}/{sample}.kegg.mapper.txt", sample = GENOMES_ANN),
                #expand("dbcan4/{sample}_dbcan/overview.txt", sample = GENOMES_DBCAN)
        output: touch("touch")

rule dbcan:
    input: "protein_fasta_dbcan/{sample}.proteins.fa"
    output: "dbcan4/{sample}_dbcan/overview.txt"
    shell:
                "run_dbcan {input} protein --out_dir dbcan4/{wildcards.sample}_dbcan --db_dir /data/ccallen/bin/run_dbcan_4.0.0/db/ --use_signalP True"

rule kofam_scan:
        input: "protein_fasta/{sample}.proteins.fa"
        output: 
                o1="kegg/{sample}/{sample}.kegg.mapper.txt",
                o2="kegg/{sample}/{sample}.ko.txt",
        shell:
                "/data/ccallen/bin/kofam_scan/exec_annotation -o {output.o2} {input} -f detail-tsv --tmp-dir tmp_d_{wildcards.sample};"
                "/data/ccallen/bin/kofam_scan/exec_annotation -o {output.o1} {input} -f mapper --tmp-dir tmp_m_{wildcards.sample};"
```

# to execute snakefile
```
cd /data/ccallen/2023_02_Sporothrix/10_annotations
conda activate run_dbcan_4.0.0 #to get an environment with snakemake installed
snakemake -n -r #to check
snakemake --cores 8
```
# 2024.04.24
* count the number of KEGG identifiers in each genome

```
cd ~/Documents/2023_02_Sporothrix/10_annotations
for f in $(tail -n +2 genomes.txt)
do
Rscript ../scripts/10_KEGG_count_identifiers.R "$f"
done
```

# 10_KEGG_count_identifiers.R
```
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

genome <- args[1]

## This assumes a fixed dir structure

library(tidyverse)
library(dplyr)

#genome = "Sporothrix_bragantina_S_CBS47491_GCA_XXXXXXXXX"

setwd(paste0("~/Documents/2023_02_Sporothrix/10_annotations/kegg/kegg_output/",genome, sep=""))
getwd()
kegg_mapper <- read.delim(paste(genome,".kegg.mapper.txt", sep=""), header = F, na.strings='') %>%
  select(-V1) %>%
  filter(!is.na(V2)) %>%
  arrange(V2)

kegg_mapper_count <- kegg_mapper %>% count(V2)
colnames(kegg_mapper_count) <- c(paste(genome), paste("n_",genome, sep=""))
write.table(kegg_mapper_count, paste(genome,".count_KEGG.txt", sep=""), sep='\t', quote = F, row.names = F)

```
* 2024.04.24
* added KEGG annotations to funannotate annotation tables

```
cd 10_annotations
for f in $(cat genomes.txt); do Rscript ../scripts/10_merge_annotations.R "$f"; done
```

* 2024.04.23
* get the list of modules and then the module files from the kegg API

```
cd /Users/carmenallen/Documents/2023_02_Sporothrix/10_annotations/kegg
wget https://rest.kegg.jp/list/module
cd mf
while IFS=$'\t' read -r v1 v2; do wget https://rest.kegg.jp/get/"$v1"; sleep 1; done < ../module
```
```
head module

M00001	Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate
M00002	Glycolysis, core module involving three-carbon compounds
M00003	Gluconeogenesis, oxaloacetate => fructose-6P
M00004	Pentose phosphate pathway (Pentose phosphate cycle)
M00005	PRPP biosynthesis, ribose 5P => PRPP
M00006	Pentose phosphate pathway, oxidative phase, glucose 6P => ribulose 5P
M00007	Pentose phosphate pathway, non-oxidative phase, fructose 6P => ribose 5P
M00008	Entner-Doudoroff pathway, glucose-6P => glyceraldehyde-3P + pyruvate
M00009	Citrate cycle (TCA cycle, Krebs cycle)
M00010	Citrate cycle, first carbon oxidation, oxaloacetate => 2-oxoglutarate
```

* calculate kegg module completeness and create heatmap using ggkegg `scripts/10_ggkegg.R`



