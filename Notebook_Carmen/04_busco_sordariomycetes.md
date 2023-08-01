# check assembly completeness using busco

* 2023.06.09
* in debary

```
/data/ccallen/2023_02_Sporothrix/04_assemblies
docker run -u $(id -u) -v $(pwd):/busco_wd ezlabgva/busco:v5.4.7_cv1 busco -i CBS139747/contigs.fasta -o CBS139747_sordario -m genome --cpu 8 	 -f

```