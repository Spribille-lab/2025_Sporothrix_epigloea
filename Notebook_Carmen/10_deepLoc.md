# use DeepLoc2 locally version 2.0 to determine the subcellular locations of each predicted protein
# 2023.07.03

```
cd /data/ccallen/2023_02_Sporothrix/10_annotations/
for sample in $(tail -n +2 genomes.txt)
do
deeploc2 -f /data/ccallen/2023_02_Sporothrix/10_annotations/protein_fasta/"$sample".proteins.fa -o DeepLoc2/"$sample" --model Accurate &>> DeepLoc2/"$sample".out
cp DeepLoc2/"$sample"/*.csv DeepLoc2/"$sample".deeploc2.csv
done
```

# well that didn't work

```
cd /data/ccallen/2023_02_Sporothrix/10_annotations/DeepLoc2
conda activate deeploc-2.0
deeploc2 -f ../protein_fasta/Sporothrix_bragantina_S_CBS47491_GCA_XXXXXXXXX.proteins.fa -o Sporothrix_bragantina_S_CBS47491_GCA_XXXXXXXXX --model Accurate &>> Sporothrix_bragantina_S_CBS47491_GCA_XXXXXXXXX.out
```

cd /data/ccallen/2023_02_Sporothrix/10_annotations/DeepLoc2
conda activate deeploc-2.0
for sample in "Sporothrix_eucalyptigena_S_CBS139899_GCA_XXXXXXXXX" "Sporothrix_eucalyptigena_S_CBS140593_GCA_XXXXXXXXX" "Sporothrix_thermara_S_CBS139747_GCA_XXXXXXXXX" "Leptographium_lundbergii_S_CBS138716_GCA_001455505" "Ophiostoma_ips_S_VPRI43529_GCA_019925475" "Ophiostoma_novoulmi_S_H327_GCA_000317715"
do
deeploc2 -f ../protein_fasta/"$sample".proteins.fa -o "$sample" --model Accurate &>> "$sample".out
done
