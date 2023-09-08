#!/bin/bash
#SBATCH --account=def-tspribi
#SBATCH --time=3-0:00
#SBATCH --cpus-per-task=8
#SBATCH --job-name=compare
#SBATCH --output=funannotate.script.Aug4.logs.out
#SBATCH --mem=249G
#SBATCH --mail-user=w6p9c9j6t9c6a2i6@spribillelabworkspace.slack.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

module load StdEnv/2020
module load apptainer/1.1.8

cd /home/ccallen/scratch/2023_02_Sporothrix/10_annotations/funannotate
apptainer run --bind $(pwd):/input $SCRATCH/docker_builds/funannotate_latest.sif funannotate compare -i \
/input/gbk/Leptographium_lundbergii_CBS138716.gbk \
/input/gbk/Ophiostoma_fasciatum_VPRI43845.gbk \
/input/gbk/Ophiostoma_ips_VPRI43529.gbk \
/input/gbk/Ophiostoma_novoulmi_H327.gbk \
/input/gbk/Sporothrix_bragantina_CBS47491.gbk \
/input/gbk/Sporothrix_brasiliensis_5110.gbk \
/input/gbk/Sporothrix_brunneoviolacea_CBS124561.gbk \
/input/gbk/Sporothrix_curviconia_CBS95973.gbk \
/input/gbk/Sporothrix_dimorphospora_CBS55374.gbk \
/input/gbk/Sporothrix_epigloea_CBS119000.gbk \
/input/gbk/Sporothrix_epigloea_CBS57363.gbk \
/input/gbk/Sporothrix_epigloea_TF4163.gbk \
/input/gbk/Sporothrix_eucalyptigena_CBS139899.gbk \
/input/gbk/Sporothrix_eucalyptigena_CBS140593.gbk \
/input/gbk/Sporothrix_euskadiensis_VPRI43754.gbk \
/input/gbk/Sporothrix_globosa_CBS120340.gbk \
/input/gbk/Sporothrix_humicola_CBS118129.gbk \
/input/gbk/Sporothrix_inflata_CBS23968.gbk \
/input/gbk/Sporothrix_insectorum_RCEF264.gbk \
/input/gbk/Sporothrix_luriei_CBS93772.gbk \
/input/gbk/Sporothrix_mexicana_CBS120341.gbk \
/input/gbk/Sporothrix_nigrograna_VPRI43755.gbk \
/input/gbk/Sporothrix_pallida_CBS13156.gbk \
/input/gbk/Sporothrix_phasma_CBS119721.gbk \
/input/gbk/Sporothrix_protearum_CBS116654.gbk \
/input/gbk/Sporothrix_pseudoabietina_VPRI43531.gbk \
/input/gbk/Sporothrix_schenckii_1099.gbk \
/input/gbk/Sporothrix_thermara_CBS139747.gbk \
/input/gbk/Sporothrix_variecibatus_CBS121961.gbk \
--cpus 8

pwd
ls .
ls $SLURM_TMPDIR
cp -r $SLURM_TMPDIR/funannotate_compare $SCRATCH/2023_02_Sporothrix/10_annotations/funannotate/.
cp -r $SLURM_TMPDIR/funannotate_compare.tar.gz $SCRATCH/2023_02_Sporothrix/10_annotations/funannotate/.
