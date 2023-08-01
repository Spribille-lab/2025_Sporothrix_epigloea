for genome in $(tail -n +2 ../genomes.txt)
do
bash ../scripts/10_gather_antismash.sh "$genome"
done

##############################################################################################

genome=$1 #Ophiostoma_fasciatum_S_VPRI43845_GCA_019925495
SP=${genome%%_S_*} #Ophiostoma_fasciatum
SPECIES=${SP//_/ } #Ophiostoma fasciatum
G=${SP%%_*} #Ophiostoma
GS=${G:0:1} #O
S=${SP#*_} #fasciatum
ss=${S:0:2} #fa
SS=$(tr '[a-z]' '[A-Z]' <<< $ss) #FA
SP_STRAIN=${genome%%_G*} #Ophiostoma_fasciatum_S_VPRI43845
STRAIN=${SP_STRAIN##*_} #VPRI43845
BASENAME="$GS""$SS""$STRAIN" #TMEATCC28783

cd /data/ccallen/2023_02_Sporothrix/10_annotations/funannotate
echo copying antiSMASH index html from "$genome"
cp "$genome"/"$BASENAME"_preds/annotate_misc//index.html antiSMASH/"$BASENAME".index.html