
KinFin


## run InterProScan on all proteins included in orthology inference

* same database as funannotate InterProScan-5.63-95.0
```
project_dir=/data/ccallen/2023_02_Sporothrix

cd "$project_dir"/09_orthology/orthofinder_input_core
cat *.fa > all_Sporothrix_proteins.fa

#run InterProScan
conda activate funannotate-1.8.15-env #to get java to work
cd "$project_dir"/09_orthology/KinFin
/data/databases/interproscan-5.63-95.0/interproscan.sh -i "$project_dir"/09_orthology/orthofinder_input_core/all_Sporothrix_proteins.fa -d out/ -dp -t p --goterms -appl Pfam -f TSV
conda deactivate

#covert to table readable by KinFin
conda activate kinfin
/data/ccallen/bin/kinfin/scripts/iprs2table.py -i out/all_Sporothrix_proteins.fa.tsv
```

0: Sporothrix_bragantina_CBS47491.fa
1: Sporothrix_brasiliensis_5110.fa
2: Sporothrix_curviconia_CBS95973.fa
3: Sporothrix_dimorphospora_CBS55374.fa
4: Sporothrix_epigloea_CBS119000.fa
5: Sporothrix_epigloea_CBS57363.fa
6: Sporothrix_epigloea_TF4163.fa
7: Sporothrix_eucalyptigena_CBS139899.fa
8: Sporothrix_eucalyptigena_CBS140593.fa
9: Sporothrix_euskadiensis_VPRI43754.fa
10: Sporothrix_globosa_CBS120340.fa
11: Sporothrix_humicola_CBS118129.fa
12: Sporothrix_inflata_CBS23968.fa
13: Sporothrix_luriei_CBS93772.fa
14: Sporothrix_mexicana_CBS120341.fa
15: Sporothrix_nigrograna_VPRI43755.fa
16: Sporothrix_pallida_CBS13156.fa
17: Sporothrix_phasma_CBS119721.fa
18: Sporothrix_protearum_CBS116654.fa
19: Sporothrix_pseudoabietina_VPRI43531.fa
20: Sporothrix_schenckii_1099.fa
21: Sporothrix_thermara_CBS139747MAG.fa
22: Sporothrix_variecibatus_CBS121961.fa


## Run KinFin
* v1.1.1
* Laetsch & Blaxter 2017

```
project_dir=/data/ccallen/2023_02_Sporothrix
cd "$project_dir"/09_orthology/KinFin
cp "$project_dir"/09_orthology/Results_Aug07/WorkingDirectory/SequenceIDs.txt .
cp "$project_dir"/09_orthology/Results_Aug07/WorkingDirectory/SpeciesIDs.txt .
cp "$project_dir"/09_orthology/Results_Aug07/Orthogroups/Orthogroups.txt .

cd "$project_dir"/09_orthology/KinFin
echo '#IDX,TAXON' > config.txt
sed 's/: /,/g' SpeciesIDs.txt | cut -f 1 -d"." >> config.txt
```
#IDX,TAXON
0,Sporothrix_bragantina_CBS47491
1,Sporothrix_brasiliensis_5110
2,Sporothrix_curviconia_CBS95973
3,Sporothrix_dimorphospora_CBS55374
4,Sporothrix_epigloea_CBS119000
5,Sporothrix_epigloea_CBS57363
6,Sporothrix_epigloea_TF4163
7,Sporothrix_eucalyptigena_CBS139899
8,Sporothrix_eucalyptigena_CBS140593
9,Sporothrix_euskadiensis_VPRI43754
10,Sporothrix_globosa_CBS120340
11,Sporothrix_humicola_CBS118129
12,Sporothrix_inflata_CBS23968
13,Sporothrix_luriei_CBS93772
14,Sporothrix_mexicana_CBS120341
15,Sporothrix_nigrograna_VPRI43755
16,Sporothrix_pallida_CBS13156
17,Sporothrix_phasma_CBS119721
18,Sporothrix_protearum_CBS116654
19,Sporothrix_pseudoabietina_VPRI43531
20,Sporothrix_schenckii_1099
21,Sporothrix_thermara_CBS139747MAG
22,Sporothrix_variecibatus_CBS121961

#IDX,TAXON,compare
0,Sporothrix_bragantina_CBS47491,other
1,Sporothrix_brasiliensis_5110,other
2,Sporothrix_curviconia_CBS95973,other
3,Sporothrix_dimorphospora_CBS55374,other
4,Sporothrix_epigloea_CBS119000,epigloea
5,Sporothrix_epigloea_CBS57363,epigloea
6,Sporothrix_epigloea_TF4163,epigloea
7,Sporothrix_eucalyptigena_CBS139899,other
8,Sporothrix_eucalyptigena_CBS140593,other
9,Sporothrix_euskadiensis_VPRI43754,other
10,Sporothrix_globosa_CBS120340,other
11,Sporothrix_humicola_CBS118129,other
12,Sporothrix_inflata_CBS23968,other
13,Sporothrix_luriei_CBS93772,other
14,Sporothrix_mexicana_CBS120341,other
15,Sporothrix_nigrograna_VPRI43755,other
16,Sporothrix_pallida_CBS13156,other
17,Sporothrix_phasma_CBS119721,other
18,Sporothrix_protearum_CBS116654,other
19,Sporothrix_pseudoabietina_VPRI43531,other
20,Sporothrix_schenckii_1099,other
21,Sporothrix_thermara_CBS139747MAG,other
22,Sporothrix_variecibatus_CBS121961,other


```
conda activate kinfin
cd "$project_dir"/09_orthology/KinFin
/data/ccallen/bin/kinfin/kinfin -g Orthogroups.txt -c config.txt -s SequenceIDs.txt -p SpeciesIDs.txt -f functional_annotation.txt -a protein_fasta/
```


network representation of clustering
```
cd "$project_dir"/09_orthology/KinFin
/data/ccallen/bin/kinfin/scripts/generate_network.py -m kinfin_results/TAXON/TAXON.cluster_summary.txt -c config.txt --exclude_universal

```



infer functional annotation of clusters based on the proteins they contain

```
cd "$project_dir"/09_orthology/KinFin
/data/ccallen/bin/kinfin/scripts/functional_annotation_of_clusters.py all -f kinfin_results/cluster_domain_annotation.IPR.txt -c kinfin_results/cluster_counts_by_taxon.txt --domain_protein_cov 0.3 --domain_taxon_cov 0.02 -o IPR
/data/ccallen/bin/kinfin/scripts/functional_annotation_of_clusters.py all -f kinfin_results/cluster_domain_annotation.GO.txt -c kinfin_results/cluster_counts_by_taxon.txt --domain_protein_cov 0.3 --domain_taxon_cov 0.02 -o GO
/data/ccallen/bin/kinfin/scripts/functional_annotation_of_clusters.py all -f kinfin_results/cluster_domain_annotation.Pfam.txt -c kinfin_results/cluster_counts_by_taxon.txt --domain_protein_cov 0.3 --domain_taxon_cov 0.02 -o Pfam
```

inferring protein count by proteome for a list of cluster IDs or protein/gene IDs

```
/data/ccallen/bin/kinfin/scripts/get_count_matrix.py --groups Orthogroups.txt --config config.txt --sequence_ids SequenceIDs.txt --functional_annotation functional_annotation.txt --cluster_ids 
