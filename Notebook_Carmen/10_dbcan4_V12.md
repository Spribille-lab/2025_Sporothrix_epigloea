# Use dbcan4.0.0 with signalP 4.1 to find CAZymes in fungal genomes and MAGs
* version 12 databases
* This tool uses three different methods of finding CAZymes
* [github for run_dbcan](https://github.com/linnabrown/run_dbcan)
* Han Zhang, Tanner Yohe, Le Huang, Sarah Entwistle, Peizhi Wu, Zhenglu Yang, Peter K Busk, Ying Xu, Yanbin Yin; dbCAN2: a meta server for automated carbohydrate-active enzyme annotation, Nucleic Acids Research, Volume 46, Issue W1, 2 July 2018, Pages W95–W101, https://doi.org/10.1093/nar/gky418
* Jinfang Zheng, Qiwei Ge, Yuchen Yan, Xinpeng Zhang, Le Huang, Yanbin Yin, dbCAN3: automated carbohydrate-active enzyme and substrate annotation, Nucleic Acids Research, 2023;, gkad328, https://doi.org/10.1093/nar/gkad328
* 2023.10.03


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
    && wget https://bcb.unl.edu/dbCAN2/download/Databases/V12/CAZyDB.07262023.fa && diamond makedb --in CAZyDB.07262023.fa -d CAZy \
    && wget https://bcb.unl.edu/dbCAN2/download/Databases/V12/dbCAN-HMMdb-V12.txt && mv dbCAN-HMMdb-V12.txt dbCAN.txt && hmmpress dbCAN.txt \
    && wget https://bcb.unl.edu/dbCAN2/download/Databases/V11/tcdb.fa && diamond makedb --in tcdb.fa -d tcdb \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/V11/tf-1.hmm && hmmpress tf-1.hmm \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/V11/tf-2.hmm && hmmpress tf-2.hmm \
    && wget https://bcb.unl.edu/dbCAN2/download/Databases/V11/stp.hmm && hmmpress stp.hmm \
    && cd ../ && wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.fna \
    && wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.faa \
    && wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.gff
```

* copy the predicted protein fastas from funannotate 1.8.5
```
cd /data/ccallen/2023_02_Sporothrix/10_annotations
for sample in $(tail -n +2 genomes_dbcan.txt) #skip the header
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
* some of the gene protein prediction was done in 2022_09_AA14/01_funannotate

```cat genomes_dbcan.txt
genome
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
Peniophora_rufa_S_SEMC34_GCA_944472145
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
Tremella_yokohamensis_S_SEMC49_GCA_944472115
Tremella_fuciformis_S_tr26_GCA_000987905
Tremella_mesenterica_S_ATCC28783_GCA_004117975
Tremella_mesenterica_S_DSM1558_GCA_000271645
Tremella_mesenterica_S_SEMC45_GCA_944472125
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

genome_df = pd.read_table('genomes_dbcan.txt').set_index("genome", drop=False)
GENOMES = list(genome_df['genome'])

rule all:
        input:
                #expand("{sample}/{sample}.kegg.mapper.txt", sample = MAGS),
                expand("dbcan4/{sample}_dbcan/overview.txt", sample = GENOMES)
                
        output: touch("touch")

rule dbcan:
    input: "protein_fasta_dbcan/{sample}.proteins.fa"
    output: "dbcan4/{sample}_dbcan/overview.txt"
    shell:
                "run_dbcan {input} protein --out_dir dbcan4/{wildcards.sample}_dbcan --db_dir /data/ccallen/bin/run_dbcan_4.0.0/db/ --use_signalP True"
```

Processed the results using R.  Used only the hmmer results 

2023_02_Sporothrix/scripts/10_process_dbcan.R

