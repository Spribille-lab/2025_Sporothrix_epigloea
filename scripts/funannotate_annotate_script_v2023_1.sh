#!/bin/bash
###################################################################################################
#                                            Funannotate Script


#This script is set up to run one assembly from start to finish.                                  #
#this script assumes that you have:
	#template.sbt.txt in the project folder
# Search for the word 'Change' to see what you will probably need to change for your set-up
##############################################################################################

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

# cp -r "$DIRNAME"/"$BASENAME"_preds ../backup #Good idea to backup somewhere
# mv "$DIRNAME"/"$BASENAME"_preds ../ANNOTATED #David moves his prediction file to a separate location.
###################################################################################################