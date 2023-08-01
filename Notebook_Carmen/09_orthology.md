



# Run fungal genomes through funannotate version 1.8.5 to get predicted proteins
* when 2023.06.16
* where /data/ccallen/2022_04_eukaryote_annotation

* how funannotate version 1.8.5 via debary fun-env

```
conda activate fun-env
export FUNANNOTATE_DB=/data/databases/funannotate_db
export PATH=$PATH:/home/ccallen/.local/bin/signalp-4.1
export GENEMARK_PATH=/home/ccallen/.local/bin/gmes_linux_64

cd /data/ccallen/2022_04_fungal_annotation/01_funannotate
bash ../scripts/funannotate_predict_script_v2022_2.sh Annulohypoxylon_stygium_S_MG137_GCA_003049155.fna &>> Annulohypoxylon_stygium_S_MG137_GCA_003049155.out
bash ../scripts/funannotate_script_v2021_1.sh Annulohypoxylon_truncatum_S_CBS140778_GCA_902805465.fna &>> Annulohypoxylon_truncatum_S_CBS140778_GCA_902805465.out

```
# for loop to loop through all the genomes.  I started using phobius locally due to pervasive problems with the remote version.
```
conda activate fun-env
export FUNANNOTATE_DB=/data/databases/funannotate_db
export PATH=$PATH:/home/ccallen/.local/bin/signalp-4.1
export GENEMARK_PATH=/home/ccallen/.local/bin/gmes_linux_64

cd /data/ccallen/2022_04_fungal_annotation/01_funannotate
for sample in $(cat samples.txt)
do
bash ../scripts/funannotate_predict_script_v2022_2.sh "$sample".fna &>> "$sample".out
done
```

# Orthology analysis of Sporothrix genomes


## Install orthofinder in debary

* 2023.06.22
```
conda create --name orthofinder2.5.5
conda activate orthofinder2.5.5
conda install orthofinder

orthofinder -h
OrthoFinder version 2.5.5 Copyright (C) 2014 David Emms

SIMPLE USAGE:
Run full OrthoFinder analysis on FASTA format proteomes in <dir>
  orthofinder [options] -f <dir>

Add new species in <dir1> to previous run in <dir2> and run new analysis
  orthofinder [options] -f <dir1> -b <dir2>

OPTIONS:
 -t <int>        Number of parallel sequence search threads [Default = 32]
 -a <int>        Number of parallel analysis threads
 -d              Input is DNA sequences
 -M <txt>        Method for gene tree inference. Options 'dendroblast' & 'msa'
                 [Default = dendroblast]
 -S <txt>        Sequence search program [Default = diamond]
                 Options: blast, diamond, diamond_ultra_sens, blast_gz, mmseqs, blast_nucl
 -A <txt>        MSA program, requires '-M msa' [Default = mafft]
                 Options: mafft, muscle
 -T <txt>        Tree inference method, requires '-M msa' [Default = fasttree]
                 Options: fasttree, raxml, raxml-ng, iqtree
 -s <file>       User-specified rooted species tree
 -I <int>        MCL inflation parameter [Default = 1.5]
 --fewer-files   Only create one orthologs file per species
 -x <file>       Info for outputting results in OrthoXML format
 -p <dir>        Write the temporary pickle files to <dir>
 -1              Only perform one-way sequence search
 -X              Don't add species names to sequence IDs
 -y              Split paralogous clades below root of a HOG into separate HOGs
 -z              Don't trim MSAs (columns>=90% gap, min. alignment length 500)
 -n <txt>        Name to append to the results directory
 -o <txt>        Non-default results directory
 -h              Print this help text

WORKFLOW STOPPING OPTIONS:
 -op             Stop after preparing input files for BLAST
 -og             Stop after inferring orthogroups
 -os             Stop after writing sequence files for orthogroups
                 (requires '-M msa')
 -oa             Stop after inferring alignments for orthogroups
                 (requires '-M msa')
 -ot             Stop after inferring gene trees for orthogroups 

WORKFLOW RESTART COMMANDS:
 -b  <dir>         Start OrthoFinder from pre-computed BLAST results in <dir>
 -fg <dir>         Start OrthoFinder from pre-computed orthogroups in <dir>
 -ft <dir>         Start OrthoFinder from pre-computed gene trees in <dir>

LICENSE:
 Distributed under the GNU General Public License (GPLv3). See License.md

CITATION:
 When publishing work that uses OrthoFinder please cite:
 Emms D.M. & Kelly S. (2019), Genome Biology 20:238

 If you use the species tree in your work then please also cite:
 Emms D.M. & Kelly S. (2017), MBE 34(12): 3267-3278
 Emms D.M. & Kelly S. (2018), bioRxiv https://doi.org/10.1101/267914
```


## put all protein fasta files into /data/ccallen/2023_02_Sporothrix/09_orthology
* 2023.07.02

input files that I shortened:
Sporothrix_bragantina_47491.proteins.fa
Sporothrix_brasiliensis_S_5110_GCA_000820605.proteins.fa
Sporothrix_brunneoviolacea_S_CBS124561_GCA_021396205.proteins.fa
Sporothrix_dimorphospora_S_CBS55374_GCA_021397985.proteins.fa
Sporothrix_epigloea_S_CBS119000_GCA_943908295.proteins.fa
Sporothrix_epigloea_S_CBS57363_GCA_943900835.proteins.fa
Sporothrix_epigloea_S_TF4163_GCA_944036445.proteins.fa
Sporothrix_eucalyptigena_CBS139899.proteins.fa
Sporothrix_eucalyptigena_CBS140593.proteins.fa
Sporothrix_euskadiensis_S_VPRI43754_GCA_019925375.proteins.fa
Sporothrix_globosa_S_CBS120340_GCA_001630435.proteins.fa
Sporothrix_humicola_S_CBS118129_GCA_021396245.proteins.fa
Sporothrix_inflata_S_CBS23968_GCA_021396225.proteins.fa
Sporothrix_insectorum_S_RCEF264_GCA_001636815.proteins.fa
Sporothrix_luriei_S_CBS93772_GCA_021398005.proteins.fa
Sporothrix_mexicana_S_CBS120341_GCA_021396375.proteins.fa
Sporothrix_nigrograna_S_VPRI43755_GCA_019925305.proteins.fa
Sporothrix_pallida_S_CBS13156_GCA_021396235.proteins.fa
Sporothrix_phasma_S_CBS119721_GCA_011037845.proteins.fa
Sporothrix_protearum_S_CBS116654_GCA_016097115.proteins.fa
Sporothrix_pseudoabietina_S_VPRI43531_GCA_019925295.proteins.fa
Sporothrix_schenckii_S_1099_GCA_000961545.proteins.fa
Sporothrix_thermara_CBS139747.proteins.fa #CBS139747_1_2_3_21_23_26_34
Sporothrix_variecibatus_S_CBS121961_GCA_021396255.proteins.fa

```
cd /data/ccallen/2023_02_Sporothrix/09_orthology
orthofinder -t 28 -f orthofinder_input
```

* also ran orthofinder using apptainer in narval so that I can try the msa and iqtree options
```
#!/bin/bash
#SBATCH --account=def-tspribi
#SBATCH --time=0-1:00
#SBATCH --cpus-per-task=8
#SBATCH --job-name=ortho_test
#SBATCH --output=orthofinder.script.logs.out
#SBATCH --mem=32G
#SBATCH --mail-user=w6p9c9j6t9c6a2i6@spribillelabworkspace.slack.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

module load StdEnv/2020
module load apptainer/1.1.8

cd $SLURM_TMPDIR
cp -r /home/ccallen/scratch/2023_02_Sporothrix/09_orthology/* .
ls .
#apptainer run --bind $(pwd):/input $SCRATCH/orthofinder_latest.sif orthofinder -o /input/orthofinder_output -f /input/orthofinder_input 

apptainer run --bind $(pwd):/input $SCRATCH/orthofinder_latest.sif orthofinder -M msa -A mafft -T iqtree -o /input/orthofinder_output -f /input/orthofinder_input
```






