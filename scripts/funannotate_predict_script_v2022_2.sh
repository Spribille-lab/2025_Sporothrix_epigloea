#!/bin/bash
###################################################################################################
#                                            Funannotate Script


#This script is set up to run one assembly from start to finish.                                  #
#this script assumes that you have:
	#Projectfolder/00_RawData for the fna files
	#projectfolder/01_funannotate for the funannotate pipeline
	#template."$BASENAME".sbt.txt in the projectfolder
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

eval "$(conda shell.bash hook)"
conda deactivate
eval "$(conda shell.bash hook)"
conda activate fun-env #version 1.8.5.  Change environment name if needed.
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
#funannotate iprscan -i "$DIRNAME"/"$BASENAME"_preds -m local --iprscan_path /data/databases/interproscan/interproscan-5.48-83.0/interproscan.sh -c 12
#funannotate remote -i "$DIRNAME"/"$BASENAME"_preds -m antismash -e ccallen@ualberta.ca & #Change to your email address. Put to the background so emapper can run at the same time.
#P1=$! #$! inherits the process ID of the previous step
#eval "$(conda shell.bash hook)"
#conda deactivate
#eval "$(conda shell.bash hook)"
#conda activate eggnog-mapper-1.0.3-env #Change to your own emapper envionment
#emapper.py  -i "$DIRNAME"/"$BASENAME"_preds/predict_results/"$SP".proteins.fa --output "$DIRNAME"/"$BASENAME"_preds/eggnog_results -d euk --data_dir /data/databases/eggnogdb --cpu 12 --override -m diamond &&
#wait $P1 #wait until process ID stops running
#eval "$(conda shell.bash hook)"
#conda deactivate
#eval "$(conda shell.bash hook)"
#conda activate fun-env
#export PATH=$PATH:/data/ccallen/bin/phobius #Change to your path
#phobius.pl -short "$DIRNAME"/"$BASENAME"_preds/predict_results/"$SP".proteins.fa > "$DIRNAME"/"$BASENAME"_preds/annotate_misc/phobius.results.txt
#funannotate annotate -i "$DIRNAME"/"$BASENAME"_preds --sbt ../template."$BASENAME".sbt.txt --eggnog "$DIRNAME"/"$BASENAME"_preds/eggnog_results.emapper.annotations --busco_db dikarya --strain "$STRAIN" --cpus 12 #Change location of template file if needed.  You may want to choose a bifferent --busco_db. You may not need the --strain flag.
# cp -r "$DIRNAME"/"$BASENAME"_preds ../backup #Good idea to backup somewhere
# mv "$DIRNAME"/"$BASENAME"_preds ../ANNOTATED #David moves his prediction file to a separate location.
###################################################################################################