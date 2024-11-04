#!/bin/bash
#SBATCH --account=def-xxxxxxx
#SBATCH --time=4-0:00
#SBATCH --cpus-per-task=32
#SBATCH --job-name=compare
#SBATCH --output=funannotate.script.Aug9.notree.logs.out
#SBATCH --mem=249G

module load StdEnv/2020
module load apptainer/1.1.8

cd $SLURM_TMPDIR
mkdir -p run_compare
cd run_compare
cp -r /scratch/ccallen/2023_02_Sporothrix/10_annotations/funannotate/gbk .

apptainer run --bind $(pwd) $SCRATCH/docker_builds/funannotate_latest.sif funannotate compare -i \
gbk/Leptographium_lundbergii_CBS138716.gbk \
gbk/Ophiostoma_fasciatum_VPRI43845.gbk \
gbk/Ophiostoma_ips_VPRI43529.gbk \
gbk/Ophiostoma_novoulmi_H327.gbk \
gbk/Sporothrix_bragantina_CBS47491.gbk \
gbk/Sporothrix_brasiliensis_5110.gbk \
gbk/Sporothrix_brunneoviolacea_CBS124561.gbk \
gbk/Sporothrix_curviconia_CBS95973.gbk \
gbk/Sporothrix_dimorphospora_CBS55374.gbk \
gbk/Sporothrix_epigloea_CBS119000.gbk \
gbk/Sporothrix_epigloea_CBS57363.gbk \
gbk/Sporothrix_epigloea_TF4163.gbk \
gbk/Sporothrix_eucalyptigena_CBS139899.gbk \
gbk/Sporothrix_eucalyptigena_CBS140593.gbk \
gbk/Sporothrix_euskadiensis_VPRI43754.gbk \
gbk/Sporothrix_globosa_CBS120340.gbk \
gbk/Sporothrix_humicola_CBS118129.gbk \
gbk/Sporothrix_inflata_CBS23968.gbk \
gbk/Sporothrix_insectorum_RCEF264.gbk \
gbk/Sporothrix_luriei_CBS93772.gbk \
gbk/Sporothrix_mexicana_CBS120341.gbk \
gbk/Sporothrix_nigrograna_VPRI43755.gbk \
gbk/Sporothrix_pallida_CBS13156.gbk \
gbk/Sporothrix_phasma_CBS119721.gbk \
gbk/Sporothrix_protearum_CBS116654.gbk \
gbk/Sporothrix_pseudoabietina_VPRI43531.gbk \
gbk/Sporothrix_schenckii_1099.gbk \
gbk/Sporothrix_thermara_CBS139747.gbk \
gbk/Sporothrix_variecibatus_CBS121961.gbk \
--num_orthos 1 --cpus 32 -o funannotate_compare_notree

cp *.tar.gz $SCRATCH/2022_03_Sporothrix/10_annotations/funannotate/.
cp funannotate_compare_notree $SCRATCH/2023_02_Sporothrix/10_annotations/funannotate/.