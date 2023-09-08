# download genome assemblies of Sporothrix using [metagenomics wiki on genome download](http://www.metagenomics.wiki/tools/fastq/ncbi-ftp-genome-download)

* 2023.08.02
* /data/ccallen/2023_02_Sporothrix


* NCBI ftp genome download
* Download list of all available genomes

```
# download list of all available genomes (GenBank), may include bad quality genomes
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt

```
* Search for available genomes of a species
* Example: Eubacterium rectale  (RefSeq database, check columns 8,9,14,15,16)

```

grep -E 'Sporothrix' assembly_summary_genbank.txt | cut -f 8,9,14,15,16

```

* Get FTP download link
* for selected genomes (e.g. Eubacterium rectale), get NCBI ftp download folder (column 20)
```
grep -E 'Sporothrix' assembly_summary_genbank.txt | cut -f 20 > ftp_folder.txt
head ftp_folder.txt

```
* extend download folder: create an exact genome (fna or gff) download link
```
awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget "ftpdir,file}' ftp_folder.txt > download_fna_files.sh
head download_fna_files.sh

```
* run download
```
source download_fna_files.sh
ls
```

* decompress genome files
```
gzip -d *.gz
ls
```

* get description (top line) of genome .fna files  (more metadata are in file assembly_summary_genbank.txt)
```
head -1 *.fna
```
* renamed all the files as Genus_species_S_strain_GCA_00XXXXXXX.fna
* convert masked .fasta file to all uppercase bases

# made a list of the genome names (without .fna)
```
for file in $(cat genomes.txt)
do
mv "$file".fna "$file".fasta
done
```

```
for file in $(cat genomes.txt); do awk 'BEGIN{FS="\n"}{if(!/>/){print toupper($0)}else{print $1}}' "$file".fasta > "$file".fna; done
rm *.fasta
```

