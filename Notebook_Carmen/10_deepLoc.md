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