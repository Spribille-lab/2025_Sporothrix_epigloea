(ENV) [ccallen@narval4 scratch]$ apptainer pull docker://davidemms/orthofinder
(ENV) [ccallen@narval4 scratch]$ apptainer run orthofinder_latest.sif orthofinder -h
(ENV) [ccallen@narval4 scratch]$ apptainer run --bind /home/ccallen/scratch/2023_02_Sporothrix/09_orthology:/input orthofinder_latest.sif orthofinder -f /input/orthofinder_input

move to a temporal drive 



apptainer pull orthofinder_latest.sif docker://ccallen/orthofinder:latest



apptainer pull docker://nextgenusfs/funannotate

apptainer pull funannotate_latest.sif docker://nextgenusfs/funannotate:latest