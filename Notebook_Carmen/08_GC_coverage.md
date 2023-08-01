# make GC coverage plots


* copy data from debary to local machine

```
cd /Users/carmenallen/Documents/2023_02_Sporothrix/08_GC_coverage
for sample in $(cat ../samples.txt)
do
mkdir -p "$sample"/data/
mkdir -p "$sample"/reports/
scp ccallen@142.244.110.136:/data/ccallen/2023_02_Sporothrix/06_binning/"$sample"/concoct_output/clustering_gt1000_merged.csv "$sample"/data/.
scp ccallen@142.244.110.136:/data/ccallen/2023_02_Sporothrix/06_binning/"$sample"/coverage_table.tsv "$sample"/data/.
scp ccallen@142.244.110.136:/data/ccallen/2023_02_Sporothrix/04_assemblies/"$sample"/contigs.fasta "$sample"/data/.
echo done "$sample"
done
```