# run eukcc2 on each file in the folder
```
cd /data/ccallen/2023_02_Sporothrix/10_annotations
conda activate eukcc-2.1.0
for file in /data/ccallen/2023_02_Sporothrix/10_annotations/protein_fasta/*.fa
do
fasta=${file##*/}
genome=${fasta%.proteins.fa}
mkdir eukcc2/"$genome"
eukcc single --out eukcc2/"$genome" --threads 8 protein_fasta/"$fasta" --db /data/databases/eukcc2_db/eukcc2_db_ver_1.1
done
```

```
for fasta in "Leptographium_lundbergii_S_CBS138716_GCA_001455505.proteins.fa" "Ophiostoma_ips_S_VPRI43529_GCA_019925475.proteins.fa" "Ophiostoma_novoulmi_S_H327_GCA_000317715.proteins.fa"
do
genome=${fasta%.proteins.fa}
mkdir eukcc2/"$genome"
eukcc single --out eukcc2/"$genome" --threads 8 protein_fasta/"$fasta" --db /data/databases/eukcc2_db/eukcc2_db_ver_1.1
done
```


# make eukcc2 report
```
for filename in */eukcc.csv
do
sed -n -e '/^\/data/p' $filename >> eukcc2_report.txt
done
```

# run busco in debary on each of the genomes

```
cd /data/ccallen/2023_02_Sporothrix/10_annotations/busco
for fasta in $(ls assemblies/*.fna)
do file=${fasta##*/}
genome=${file%.*}
echo $genome
docker run -u $(id -u) -v $(pwd):/busco_wd ezlabgva/busco:v5.4.7_cv1 busco -i "$fasta" -o "$genome" -m genome -l sordariomycetes --cpu 8 -f
done
```

compile report
```
cd /data/ccallen/2023_02_Sporothrix/10_annotations/busco
for filename in */*.txt
do
echo -n  $filename >> busco_report.txt;
sed -n -e '/^\tC:/p' $filename >> busco_report.txt
done
```