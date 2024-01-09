# Use dbcan4.0.0 with signalP 4.1 to find CAZymes in fungal genomes and MAGs
* version 11 databases
* This tool uses three different methods of finding CAZymes
* [github for run_dbcan](https://github.com/linnabrown/run_dbcan)
* Han Zhang, Tanner Yohe, Le Huang, Sarah Entwistle, Peizhi Wu, Zhenglu Yang, Peter K Busk, Ying Xu, Yanbin Yin; dbCAN2: a meta server for automated carbohydrate-active enzyme annotation, Nucleic Acids Research, Volume 46, Issue W1, 2 July 2018, Pages W95–W101, https://doi.org/10.1093/nar/gky418
* Jinfang Zheng, Qiwei Ge, Yuchen Yan, Xinpeng Zhang, Le Huang, Yanbin Yin, dbCAN3: automated carbohydrate-active enzyme and substrate annotation, Nucleic Acids Research, 2023;, gkad328, https://doi.org/10.1093/nar/gkad328
* 2023.07.03


05/01/2023: dbCAN3 paper is published at Nucleic Acids Research featuring substrate prediction
02/11/2023: dbCAN2 is updated to dbCAN3 with glycan substrate prediction functions: 1. CAZyme substrate prediction based on dbCAN-sub ; 2. CGC substrate prediction based on dbCAN-PUL searching and dbCAN-sub majority voting. For CGC substrate prediction, please see our dbCAN-seq update paper for details. With these new functions (esp. the dbCAN-sub search), dbCAN3 is now slower to get the result back to you. Please be patience!
8/9/2022: dbCAN HMMdb v11 is released (based on CAZyDB 8/7/2022). Now the HMMdb contains 699 CAZyme HMMs (452 family HMMs + 3 cellulosome HMMs + 244 subfamily HMMs). The CAZyDB for Diamond search is also updated, containing in total 2,428,817 fasta sequences. See readme for details.
06/29/2022: dbCAN-sub (HMMdb from eCAMI subfams and allows EC and substrate inferences) is now deployed on dbCAN meta server and replaces eCAMI (consumes too much RAM and too slow).

```

conda create --name run_dbcan_4.0.0 python=3.8 dbcan -c conda-forge -c bioconda
conda install -c conda-forge mamba
mamba install -c conda-forge -c bioconda snakemake


cd /data/ccallen/bin/run_dbcan_4.0.0
test -d db || mkdir db
cd db \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/fam-substrate-mapping-08252022.tsv \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/PUL.faa && makeblastdb -in PUL.faa -dbtype prot \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL_07-01-2022.xlsx \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL_07-01-2022.txt \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL.tar.gz && tar xvf dbCAN-PUL.tar.gz \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN_sub.hmm && hmmpress dbCAN_sub.hmm \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/V11/CAZyDB.08062022.fa && diamond makedb --in CAZyDB.08062022.fa -d CAZy \
    && wget https://bcb.unl.edu/dbCAN2/download/Databases/V11/dbCAN-HMMdb-V11.txt && mv dbCAN-HMMdb-V11.txt dbCAN.txt && hmmpress dbCAN.txt \
    && wget https://bcb.unl.edu/dbCAN2/download/Databases/V11/tcdb.fa && diamond makedb --in tcdb.fa -d tcdb \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/V11/tf-1.hmm && hmmpress tf-1.hmm \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/V11/tf-2.hmm && hmmpress tf-2.hmm \
    && wget https://bcb.unl.edu/dbCAN2/download/Databases/V11/stp.hmm && hmmpress stp.hmm \
    && cd ../ && wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.fna \
    && wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.faa \
    && wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.gff
```

* copy the predicted protein fastas from funannotate
```
cd /data/ccallen/2023_02_Sporothrix/10_annotations
for sample in $(tail -n +2 genomes.txt) #skip the header
do
SP=${sample%%_S_*} #Tremella_mesenterica
SPECIES=${SP//_/ } #Tremella mesenterica
G=${SP%%_*} #Tremella
GS=${G:0:1} #T
S=${SP#*_} #mesenterica
ss=${S:0:2} #me
SS=$(tr '[a-z]' '[A-Z]' <<< $ss) #ME
SP_STRAIN=${sample%%_G*} #Tremella_mesenterica_S_ATCC28783
STRAIN=${SP_STRAIN##*_} #ATCC28783
BASENAME="$GS""$SS""$STRAIN" #TMEATCC28783

echo working on $sample
mkdir protein_fasta/"$sample"
cp /data/ccallen/2022_04_fungal_annotation/01_funannotate/"$sample"/*_preds/predict_results/"$SP".proteins.fa protein_fasta/"$sample"/.
mv protein_fasta/"$sample"/"$SP".proteins.fa protein_fasta/"$sample".proteins.fa
rm -r protein_fasta/"$sample"
done
```

```cat genomes.txt
genome
Sporothrix_bragantina_S_CBS47491_GCA_XXXXXXXXX
Sporothrix_brasiliensis_S_5110_GCA_000820605
Sporothrix_brunneoviolacea_S_CBS124561_GCA_021396205
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
```

## run dbcan using snakemake

```
cd /data/ccallen/2023_02_Sporothrix/10_annotations
export PATH=$PATH:/home/ccallen/.local/bin/signalp-4.1
conda activate run_dbcan_4.0.0
snakemake -np #to do a dry run
snakemake --cores 8
```

```snakefile
import glob
import os
import pandas as pd

genome_df = pd.read_table('genomes.txt').set_index("genome", drop=False)
GENOMES = list(genome_df['genome'])

rule all:
        input:
                #expand("{sample}/{sample}.kegg.mapper.txt", sample = MAGS),
                expand("dbcan3/{sample}_dbcan/overview.txt", sample = GENOMES)
                #expand("DeepLoc2/{sample}_dbcan/overview.txt", sample = GENOMES)
                
        output: touch("touch")

rule dbcan:
    input: "protein_fasta/{sample}.proteins.fa"
    output: "dbcan3/{sample}_dbcan/overview.txt"
    shell:
                "run_dbcan {input} protein --out_dir dbcan3/{wildcards.sample}_dbcan --db_dir /data/ccallen/bin/run_dbcan_4.0.0/db/ --use_signalP True"
```



Processed the results using R

2023_02_Sporothrix/scripts/10_process_dbcan.R