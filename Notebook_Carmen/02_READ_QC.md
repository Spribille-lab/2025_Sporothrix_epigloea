# use metaWRAP module to clean the original reads
* using metaWRAP version 1.3.2
* 2023.04.07

```
cd /data/ccallen/2023_02_Sporothrix/01_RAW_READS
gzip -d *.fastq.gz
```

* clean the data
```
conda activate metawrap-env
export PATH=$PATH:/data/ccallen/packages/metaWRAP/bin
cd /data/ccallen/2023_02_Sporothrix
for sample in $(cat samples.txt)
do
metawrap read_qc -t 28 -1 01_RAW_READS/"$sample"_1.fastq -2 01_RAW_READS/"$sample"_2.fastq -o 02_READ_QC/"$sample"
done
```



* copy QC reports to local machine

```
cd /data/ccallen/2023_02_Sporothrix/02_READ_QC
for sample in $(cat ../samples.txt)
do
mkdir qc_html/"$sample"
cp "$sample"/post-QC_report/final_pure_reads_1_fastqc.html qc_html/"$sample"/.
mv qc_html/"$sample"/final_pure_reads_1_fastqc.html qc_html/"$sample".final_pure_reads_1_fastqc.html
cp "$sample"/post-QC_report/final_pure_reads_2_fastqc.html qc_html/"$sample"/.
mv qc_html/"$sample"/final_pure_reads_2_fastqc.html qc_html/"$sample".final_pure_reads_2_fastqc.html
rm -r qc_html/"$sample"
done
```