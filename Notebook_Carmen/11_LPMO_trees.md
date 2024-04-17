


## Processed the results of dbcan4 using R.  Used only the hmmer results 

2023_02_Sporothrix/scripts/10_process_dbcan.R

Made a table of lpmo and AA14 coordinates using the above script.


* make one fasta file with AA14/LPMO sequences of all genomes
```
cd /data/ccallen/2023_02_Sporothrix/10_annotations/dbcan4
while IFS=$'\t' read -r v1 v2 v3 v4; do samtools faidx ../protein_fasta_dbcan/"$v1".proteins.fa $v2:$v3-$v4 >> summarized_outputs/allgenomes.cazy_aa14.fa; done < summarized_outputs/AA14_coordinates.txt
while IFS=$'\t' read -r v1 v2 v3 v4; do samtools faidx ../protein_fasta_dbcan/"$v1".proteins.fa $v2:$v3-$v4 >> summarized_outputs/Sordariomycetes.cazy_aa14.fa; done < summarized_outputs/AA14_Sordariomycetes_coordinates.txt
while IFS=$'\t' read -r v1 v2 v3 v4; do samtools faidx ../protein_fasta_dbcan/"$v1".proteins.fa $v2:$v3-$v4 >> summarized_outputs/allgenomes.cazy_lpmo.fa; done < summarized_outputs/lpmo_coordinates.txt
```

renamed the fasta headers to include the family in the header name (see 10_process_dbcan.R)

## tree of lpmos extracted from annotated genomes

align lpmos using mafft in geneious prime
/Users/carmenallen/Documents/2023_02_Sporothrix/10_annotations/dbcan4/allgenomes.cazy_lpmo.fa

geneious prime Geneious Prime® 2023.2.1 Build 2023-07-20 11:29 Java Version 11.0.18+10 (64 bit) Installation Id cdAwcQoBDwdisnkNCyWctA==
MAFFT v7.450 algorithm auto, Scoring matrix BLOSUM62 Gap open penalty 1.53 Offset value 0.123

exported phylip format

```
cd /Users/carmenallen/Documents/2023_02_Sporothrix/11_LPMO_trees/LPMO
conda activate trimal_1.4.1 
trimal -in LPMO_alignment.phy -out LPMO_alignment_trimmed.phy -gappyout

cd /Users/carmenallen/Documents/2023_02_Sporothrix/11_LPMO_trees/LPMO
/Users/carmenallen/Programs/iqtree-2.0-rc1-MacOSX/bin/iqtree -s /Users/carmenallen/Documents/2023_02_Sporothrix/11_LPMO_trees/LPMO/LPMO_alignment_trimmed.phy -m MFP -B 1000 -pre LPMO
```

make a tree of lpmos from annotated genomes using iqtree

## extract Sordariomycete AA14 sequences from genbank and use those to make a tree along with AA14s from Sordariomycete genomes that I annotated.

2023.10.04
* blastx KY769370 from P. coccineus
* Your search is limited to records that include: Sordariomycetes (taxid:147550); and exclude: Fusarium (taxid:5506)
* nr database
* took 1000 results
* 1e-25 or better
* downloaded fasta (130 sequences)

2023.10.04
* blastx KY769370 from P. coccineus
* Your search is limited to records that include Fusarium (taxid:5506)
* nr database
* took 10 results
* 1e-43 or better
* downloaded fasta (10 sequences)

/Users/carmenallen/Documents/2023_02_Sporothrix/11_LPMO_trees/Sordariomycetes/blastx_sordariomycetes.fa
/Users/carmenallen/Documents/2023_02_Sporothrix/11_LPMO_trees/Sordariomycetes/blastx_fusarium.fa
rename fasta headers in sublime

(^>\S+\.\d)(.+)(\[)(\S+)(\s)(\S+)(\s*)(\S*)(\])
\1_\4_\6_\8

(^>\S+\.\d)(.+)(\[)(\S+)(\s)(\S+)(\s*)(\S*)(\s*)(\S*)(\])
\1_\4_\6_\8$10

(^>\S+\.\d)(.+)(\_$)
\1\2


## AA14 tree

2023.10.05
* blastx KY769370 from P. coccineus
* nr database
* 1e-25 or better

100 Basidiomycota sequences (top e-value)

2023.10.05
* blastx KY769370 from P. coccineus
* nr database
* 1e-25 or better
* Your search is limited to records that include: Ascomycota (taxid:4890) ; and exclude: Sordariomycetes (taxid:147550)
* 209 sequences

combined fastas of Ascomycota (minus Sordario), Basidiomycota, and Sordariomycete AA14s from blastx hits with annotated genomes

```
cd /Users/carmenallen/Documents/2023_02_Sporothrix/
cat 11_LPMO_trees/AA14/Sordariomycetes/blastx_fusarium.fa 11_LPMO_trees/AA14/Sordariomycetes/blastx_sordariomycetes.fa 11_LPMO_trees/AA14/Fungi/blastx_Ascomycota_exclude_Sordario.fa 11_LPMO_trees/AA14/Fungi/blastx_basidiomycota.fa 10_annotations/dbcan4/summarized_outputs/allgenomes.cazy_aa14.renamed.fa > 11_LPMO_trees/AA14/Fungi/Fungi_AA14.fa
```
In total 508 sequences

geneious prime Geneious Prime® 2023.2.1 Build 2023-07-20 11:29 Java Version 11.0.18+10 (64 bit) Installation Id cdAwcQoBDwdisnkNCyWctA==
MAFFT v7.450 algorithm auto, Scoring matrix BLOSUM62 Gap open penalty 1.53 Offset value 0.123

exported phylip format

```
cd /Users/carmenallen/Documents/2023_02_Sporothrix/11_LPMO_trees/AA14/Fungi
conda activate trimal_1.4.1 
trimal -in Fungi_AA14_alignment.phy -out Fungi_AA14_alignment_trimmed.phy -gappyout

cd /Users/carmenallen/Documents/2023_02_Sporothrix/11_LPMO_trees/AA14/Fungi
/Users/carmenallen/Programs/iqtree-2.0-rc1-MacOSX/bin/iqtree -s /Users/carmenallen/Documents/2023_02_Sporothrix/11_LPMO_trees/AA14/Fungi/Fungi_AA14_alignment_trimmed.phy -m MFP -B 1000 -pre AA14_Fungi
```


## Sordariomycetes AA14 tree

combined fastas of Sordariomycete AA14s from annotated genomes with blastx hits

```
cd /Users/carmenallen/Documents/2023_02_Sporothrix/
cat 11_LPMO_trees/AA14/Sordariomycetes/blastx_fusarium.fa 11_LPMO_trees/AA14/Sordariomycetes/blastx_sordariomycetes.fa 10_annotations/dbcan4/summarized_outputs/Sordariomycetes.cazy_aa14.fa > 11_LPMO_trees/AA14/Sordariomycetes/Sordariomycetes_AA14.fa
```

geneious prime Geneious Prime® 2023.2.1 Build 2023-07-20 11:29 Java Version 11.0.18+10 (64 bit) Installation Id cdAwcQoBDwdisnkNCyWctA==
MAFFT v7.450 algorithm auto, Scoring matrix BLOSUM62 Gap open penalty 1.53 Offset value 0.123

exported phylip format


```
cd /Users/carmenallen/Documents/2023_02_Sporothrix/11_LPMO_trees/AA14/Sordariomycetes
conda activate trimal_1.4.1 
trimal -in Sordariomycetes_AA14_alignment.phy -out Sordariomycetes_AA14_alignment_trimmed.phy -gappyout

cd /Users/carmenallen/Documents/2023_02_Sporothrix/11_LPMO_trees/Sordariomycetes
/Users/carmenallen/Programs/iqtree-2.0-rc1-MacOSX/bin/iqtree -s /Users/carmenallen/Documents/2023_02_Sporothrix/11_LPMO_trees/Sordariomycetes/Sordariomycetes_AA14_alignment_trimmed.phy -m MFP -B 1000 -pre AA14_sordario
```


# re-do March 2024
need to use blastp this time


## AA14 tree

2024.03.07
/Users/carmenallen/Documents/2023_02_Sporothrix/11_LPMO_trees/AA14/Fungi_March2024

query: original AA14 characterized AUM86167.1
run tree in debary with iqtree2.2.3
blastp to the new experimental database nr_clustered
evalue cutoff of 1e-25
1308 blastp hits, 59 from genome annotations


* renamed fasta file 

rename fasta headers in sublime

(^>\S+\.\d)(.+)(\[)(\S+)(\s)(\S+)(\s*)(\S*)(\])
\1_\4_\6_\8

(^>\S+\.\d)(.+)(\[)(\S+)(\s)(\S+)(\s*)(\S*)(\s*)(\S*)(\])
\1_\4_\6_\8$10

(^>\S+\.\d)(.+)(\_$)
\1\2


* combined fastas from blastp hits with annotated genomes

```
cd /Users/carmenallen/Documents/2023_02_Sporothrix/
cat 11_LPMO_trees/AA14/Fungi_March2024/blastp/blastp_results.fa 10_annotations/dbcan4/summarized_outputs/allgenomes.cazy_aa14.renamed.fa > 11_LPMO_trees/AA14/Fungi_March2024/Fungi_AA14.fa
```

aligned using geneious
Geneious Prime® 2024.0.2 Build 2024-02-14 15:15 Java Version 11.0.20.1+1 (64 bit)
MAFFT v7.450 algorithm auto, Scoring matrix BLOSUM62 Gap open penalty 1.53 Offset value 0.123

* exported phylip format (relaxed)
* trim

```
cd /Users/carmenallen/Documents/2023_02_Sporothrix/11_LPMO_trees/AA14/Fungi_March2024
conda activate trimal_1.4.1 
trimal -in Fungi_AA14_alignment.phy -out Fungi_AA14_alignment_trimmed.phy -gappyout

```
* tree
```
conda activate iqtree2-env
cd /data/ccallen/2023_02_Sporothrix/11_LPMO_trees/AA14/Fungi_March2024
iqtree -s Fungi_AA14_alignment_trimmed.phy -nt AUTO -m MFP -B 1000 -pre AA14_Fungi
```

* extracted taxonomic information using taxize 0.9.100
2024.03.13

11_taxize.R

* used itol to visualize the tree
annotation text files are in:
/Users/carmenallen/Documents/2023_02_Sporothrix/11_LPMO_trees/AA14/Fungi_March2024/itol 
