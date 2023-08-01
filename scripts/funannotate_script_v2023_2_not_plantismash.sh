#!/bin/bash
###################################################################################################
#                                            Funannotate Script


#This script is set up to run one assembly from start to finish.                                  #
#this script assumes that you have:
	#Projectfolder/00_RawData for the fna files
	#projectfolder/01_funannotate for the funannotate pipeline
	#template."$BASENAME".sbt.txt in the projectfolder
# Search for the word 'Change' to see what you will probably need to change for your set-up


# Adapted from Carmen Allen's script (10-01-2022)
##############################################################################################

# Sample filename: AC_H1932_watanabae_bin.20.fa


FILENAME="${1##*/}" # remove anything before the '/'

DIRNAME="${FILENAME%.*}" # remove the file extension (.fa)

SPECIES_="${DIRNAME#*_*_}" # remove leading info ("AC_ISOLATE_")
SPECIES_="${SPECIES_%%_bin*}" # remove everything after species name ("_bin_num")
SPECIES_="${SPECIES_^}" # make first letter capitalized (only works on bsh, not zsh)

SPECIES="${SPECIES_//_/ }" # replace the underscore with space in species name

ISOLATE="${DIRNAME#*_}" # remove leading info ("AC_")
ISOLATE="${ISOLATE%%_*}" # remove trailing info after isolate info ("_species_name...")

# make new filename
GS="${SPECIES:0:1}" # extract first letter of genus name
SS="${SPECIES#* }" # extract species name from binom
SS="${SS:0:2}" # extract first two letters of species name
SS="${SS^^}" # capitalize both species letters (only works on bsh, not zsh)
BASENAME="$GS""$SS""$ISOLATE" # combine into one basename


# Set conda environment to funannotate
eval "$(conda shell.bash hook)" # shut down whatever's running (if it's running)
conda deactivate
eval "$(conda shell.bash hook)" # activate the funannotate enviroment
conda activate funannotate-env # CHANGE: change to your environment

# Setup folder to send results
mkdir "$DIRNAME"

# Run funannotate initial steps (clean, sort, mask) to setup for predictions
funannotate clean -i ../00_rawdata/"$FILENAME" -o  "$DIRNAME"/"$DIRNAME"_cleaned.fas --exhaustive &&
funannotate sort -i  "$DIRNAME"/"$DIRNAME"_cleaned.fas -o  "$DIRNAME"/"$DIRNAME"_sorted.fas &&
funannotate mask -i  "$DIRNAME"/"$DIRNAME"_sorted.fas -o  "$DIRNAME"/"$DIRNAME"_masked.fas --cpus 12

mv ./funannotate-mask.log "$DIRNAME"/funannotate-mask.log

# run funnanotate predict
funannotate predict -i "$DIRNAME"/"$DIRNAME"_masked.fas -s "$SPECIES" --isolate "$ISOLATE" -o "$DIRNAME"/"$BASENAME"_preds --optimize_augustus --repeats2evm --max_intronlen 50000 --cpus 12 --organism other --busco_seed_species chlorella --busco_db eukaryota --name "$BASENAME"

# run iprscan in funannotate
funannotate iprscan -i "$DIRNAME"/"$BASENAME"_preds -m local --iprscan_path /data/databases/interproscan/interproscan-5.48-83.0/interproscan.sh -c 12

# run phobius
/data/andrewc/scripts/phobius/phobius.pl -short "$DIRNAME"/"$BASENAME"_preds/predict_results/"$SPECIES_"_"$ISOLATE".proteins.fa > "$DIRNAME"/"$BASENAME"_preds/phobius.results.txt

# to run it remotely
#funannotate remote -i "$DIRNAME"/"$BASENAME"_preds -m phobius -e acook@ualberta.ca

# run plantismash
#eval "$(conda shell.bash hook)"
#conda deactivate
#eval "$(conda shell.bash hook)"
#conda activate plantismash-env

#/data/andrewc/scripts/antismash/plantismash-1.0/run_antismash.py --taxon plants --logfile "$DIRNAME"/"$BASENAME"_preds/logfiles/antismash.log "$DIRNAME"/"$BASENAME"_preds/predict_results/"$SPECIES_"_"$ISOLATE".gbk

# run emapper
#eval "$(conda shell.bash hook)"
#conda deactivate
#eval "$(conda shell.bash hook)"
#conda activate funannotate-env

# run eggnog mapper
emapper.py  -i "$DIRNAME"/"$BASENAME"_preds/predict_results/"$SPECIES_"_"$ISOLATE".proteins.fa --output "$DIRNAME"/"$BASENAME"_preds/eggnog_results -d euk --data_dir /data/databases/eggnogdb5 --cpu 12 --override -m diamond

# run funannotate's annotate feature!
funannotate annotate -i "$DIRNAME"/"$BASENAME"_preds --sbt ../template.sbt.txt --eggnog "$DIRNAME"/"$BASENAME"_preds/eggnog_results.emapper.annotations --iprscan "$DIRNAME"/"$BASENAME"_preds/annotate_misc/iprscan.xml --phobius "$DIRNAME"/"$BASENAME"_preds/phobius.results.txt --busco_db eukaryota --strain "$ISOLATE" --cpus 12

# --antismash "$DIRNAME"/"$SPECIES_"_"$ISOLATE"/scaffold_1.final.gbk