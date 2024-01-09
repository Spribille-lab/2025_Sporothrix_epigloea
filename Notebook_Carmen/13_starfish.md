# prepare fasta file and gff3 files

extracted the contig headers

```
cd /data/ccallen/2023_02_Sporothrix/12_starfish/assembly
grep -e ">" "$genome".contigs.fasta > "$genome".headers
```

used sublime to make a lookup table
```
find (>)(NODE)(_)(\d+)(.+)
replace \2\3\4\5\tSEPCBS119000_\2\4
```

on a mac:
```
while IFS=$'\t' read -r old new; do gsed -i "s/$old/$new/g" Sporothrix_epigloea_S_CBS119000_GCA_943908295.contigs.fasta.renamed; done < Sporothrix_epigloea_S_CBS119000_GCA_943908295.lookup
```

in compute canada
```
genome=Sporothrix_epigloea_S_CBS57363_GCA_943900835
while IFS=$'\t' read -r old new; do sed -i "s/$old/$new/g" "$genome".contigs.final.fasta; done < "$genome".lookup
```
Sporothrix_epigloea_S_CBS57363_GCA_XXXXXXXXX.scaffolds.fasta
Sporothrix_epigloea_S_CBS57363_GCA_XXXXXXXXX.scaffolds.lookup


This should work but I couldn't get it to work
```
sed -r -i 's/(NODE\_\\d+)(.*)/CBS119000_\1/' test.fasta
```

copied gff3 files from latest annotation in debary to starfish analysis
```
cd /data/ccallen/2023_02_Sporothrix/12_starfish
cp /data/ccallen/2023_02_Sporothrix/10_annotations/funannotate/"$genome"/*_preds/annotate_results/*.gff3 gff3/.
```

rename the contigs to include the genome name
```
genome=Sporothrix_epigloea_S_CBS119000_GCA_943908295
sed 's/scaffold_/SEPCBS119000_NODE/' "$genome".gff3 | grep -v '#' > "$genome".final.gff3
genome=Sporothrix_epigloea_S_CBS57363_GCA_943900835
sed 's/scaffold_/SEPCBS57363_NODE/' "$genome".gff3 | grep -v '#' > "$genome".final.gff3
genome=Sporothrix_epigloea_S_TF4163_GCA_944036445
sed 's/scaffold_/SEPTF4163_NODE/' "$genome".gff3 | grep -v '#' > "$genome".final.gff3
```

# start tutorial

navigate to the starfish data directory:
```
cd /data/ccallen/2023_02_Sporothrix/12_starfish
```
create ome2*.txt files detailing the absolute path to each genome's gff3 and assembly:
```
realpath assembly/*.final.fasta | perl -pe 's/^(.+?([^\/]+?).contigs.final.fasta)$/\2\t\1/' > ome2assembly.txt
realpath gff3/*.final.gff3 | perl -pe 's/^(.+?([^\/]+?).final.gff3)$/\2\t\1/' > ome2gff.txt
```
concatenate all gff3 files into a single file (a useful shortcut for some analyses):
```
cat gff3/*.final.gff3 > sep3.gff3
```

concatenate all assembly files and make a blastn database:
```
mkdir blastdb
cut -f2 ome2assembly.txt | xargs cat > blastdb/spor29.assemblies.fna
makeblastdb -in blastdb/spor29.assemblies.fna -out blastdb/spor29.assemblies -parse_seqids -dbtype nucl
```

calculate %GC content across all genomes (useful for visualizing elements later):
```
/data/ccallen/bin/starfish/scripts/seq-gc.sh -Nbw 1000 blastdb/spor29.assemblies.fna > spor29.assemblies.gcContent_w1000.bed
rm blastdb/spor29.assemblies.fna
```

parse the provided eggnog mapper annotations (NB the format of the output file has changed in more recent emapper versions):
had to do this in compute canada because sed doesn't work in debary
```
cut -f1,8  ann/*emapper.annotations | grep -v  '#' | grep -v -P '\t-' | perl -pe 's/\t/\tEMAP\t/' | grep -vP '\tNA' | sed 's/-T1//' > ann/sep3.gene2emap.txt
```

In order to parse the output from more recent versions of emapper, you must use the following command to retrieve the narrowest eggnog ortholog group per sequence:
compute canada
```
cut -f1,5 ann/*emapper.annotations | grep -v '#' | perl -pe 's/(^.+?)\t.+,([^,]+)$/\1\t\2/' | perl -pe 's/@/\t/' | sed 's/-T1//'> ann/sep3.gene2og.txt
```

convert to .mcl format:
```
/data/ccallen/bin/starfish/scripts/geneOG2mclFormat.pl -i ann/macph6.gene2og.txt -o ann/
```

# Gene finder module


We begin by *de novo* annotating all tyrosine recombinases (tyrs/YRs) in the provided assemblies. In practice, we can de novo annotate any gene we want, as long as we have an HMM file of a predicted domain within that gene and a multifasta of amino acid sequences of that gene (the more predicted sequences the better).

first, create a dedicated directory for good housekeeping:
```
mkdir geneFinder
```

*de novo* annotate tyrs with the provided YR HMM and amino acid queries (~10min):
```
starfish annotate -T 28 -x spor29_tyr -a ome2assembly.txt -g ome2gff.txt -p /data/ccallen/bin/starfish/database/YRsuperfams.p1-512.hmm -P /data/ccallen/bin/starfish/database/YRsuperfamRefs.faa -i tyr -o geneFinder/
```

had to remove the first line from gff3 files 

```
genome=Sporothrix_epigloea_S_CBS57363_GCA_943900835
grep -v '#' "$genome".final.gff3 | grep 'gene' > "$genome".final.final.gff3

genome=Sporothrix_epigloea_S_CBS119000_GCA_943908295
grep -v '#' "$genome".final.gff3 | grep 'gene' > "$genome".final.final.gff3

genome=Sporothrix_epigloea_S_TF4163_GCA_944036445
grep -v '#' "$genome".final.gff3 | grep 'gene' > "$genome".final.final.gff3
```
Ended up using regex in sublime to get gff3 file formatted properly.  This is how they have to look:

```
SEPCBS119000_NODE1	funannotate	gene	662	1994	.	+	.	ID=SEPCBS119000_000001;Name=SEPCBS119000_000001
SEPCBS119000_NODE1	funannotate	gene	2871	5464	.	-	.	ID=SEPCBS119000_000002;Name=SEPCBS119000_000002
SEPCBS119000_NODE1	funannotate	gene	6846	13263	.	-	.	ID=SEPCBS119000_000003;Name=SEPCBS119000_000003
SEPCBS119000_NODE1	funannotate	gene	20451	22465	.	+	.	ID=SEPCBS119000_000004;Name=SEPCBS119000_000004
SEPCBS119000_NODE1	funannotate	gene	22864	24006	.	-	.	ID=SEPCBS119000_000005;Name=SEPCBS119000_000005
SEPCBS119000_NODE1	funannotate	gene	24824	25376	.	-	.	ID=SEPCBS119000_000006;Name=SEPCBS119000_000006
SEPCBS119000_NODE1	funannotate	gene	26150	26899	.	-	.	ID=SEPCBS119000_000007;Name=SEPCBS119000_000007
SEPCBS119000_NODE1	funannotate	gene	27702	28479	.	-	.	ID=SEPCBS119000_000008;Name=SEPCBS119000_000008
SEPCBS119000_NODE1	funannotate	gene	28886	29345	.	+	.	ID=SEPCBS119000_000009;Name=SEPCBS119000_000009
```



Sporothrix_bragantina_S_CBS47491_GCA_XXXXXXXXX.out:masked repeats: 2,123,843 bp (5.58%)
Sporothrix_brasiliensis_S_5110_GCA_000820605.out:masked repeats: 2,111,403 bp (6.36%)
Sporothrix_brunneoviolacea_S_CBS124561_GCA_021396205.out:masked repeats: 2,758,461 bp (7.33%)
Sporothrix_brunneoviolacea_S_CBS_124561_GCA_021396205.out:masked repeats: 2,758,461 bp (7.33%)
Sporothrix_curviconia_S_CBS95973_GCA_016097085.out:masked repeats: 2,592,864 bp (7.39%)
Sporothrix_curviconia_S_CBS95973_GCA_XXXXXXXXX.out:masked repeats: 2,857,770 bp (7.96%)
Sporothrix_dimorphospora_S_CBS55374_GCA_021397985.out:masked repeats: 2,893,142 bp (7.44%)
Sporothrix_epigloea_S_CBS119000_GCA_943908295.out:masked repeats: 647,187 bp (2.79%)
Sporothrix_epigloea_S_CBS57363_GCA_943900835.out:masked repeats: 623,518 bp (2.71%)
Sporothrix_epigloea_S_TF4163_GCA_944036445.out:masked repeats: 575,308 bp (2.50%)
Sporothrix_eucalyptigena_S_CBS139899_GCA_XXXXXXXXX.out:masked repeats: 2,081,892 bp (5.31%)
Sporothrix_eucalyptigena_S_CBS140593_GCA_XXXXXXXXX.out:masked repeats: 2,088,815 bp (5.32%)
Sporothrix_euskadiensis_S_VPRI43754_GCA_019925375.out:masked repeats: 2,785,353 bp (7.71%)
Sporothrix_globosa_S_CBS120340_GCA_001630435.out:masked repeats: 1,916,395 bp (5.73%)
Sporothrix_humicola_S_CBS118129_GCA_021396245.out:masked repeats: 3,081,860 bp (7.57%)
Sporothrix_inflata_S_CBS23968_GCA_021396225.out:masked repeats: 3,033,645 bp (7.68%)
Sporothrix_inflata_S_CBS23968_GCA_021396225.out:masked repeats: 3,033,645 bp (7.68%)
Sporothrix_insectorum_S_RCEF264_GCA_001636815.out:masked repeats: 3,116,887 bp (8.98%)
Sporothrix_luriei_S_CBS93772_GCA_021398005.out:masked repeats: 2,312,610 bp (7.10%)
Sporothrix_mexicana_S_CBS120341_GCA_021396375.out:masked repeats: 2,864,689 bp (6.55%)
Sporothrix_nigrograna_S_VPRI43755_GCA_019925305.out:masked repeats: 2,655,076 bp (9.71%)
Sporothrix_pallida_S_CBS13156_GCA_021396235.out:masked repeats: 2,891,946 bp (7.19%)
Sporothrix_pallida_S_SPA8_GCA_000710705.out:masked repeats: 2,734,434 bp (7.23%)
Sporothrix_phasma_S_CBS119588_GCA_016097075.out:masked repeats: 1,993,897 bp (6.44%)
Sporothrix_phasma_S_CBS119721_GCA_011037845.out:masked repeats: 2,050,657 bp (6.79%)
Sporothrix_protearum_S_CBS116654_GCA_016097115.out:masked repeats: 3,265,502 bp (8.42%)
Sporothrix_pseudoabietina_S_VPRI43531_GCA_019925295.out:masked repeats: 2,630,450 bp (7.47%)
Sporothrix_schenckii_S_1099_GCA_000961545.out:masked repeats: 1,948,555 bp (6.02%)
Sporothrix_schenckii_S_ATCC58251_GCA_000474925.out:masked repeats: 1,890,787 bp (5.87%)
Sporothrix_schenckii_S_SsEM7_GCA_002837075.out:masked repeats: 1,999,560 bp (6.08%)
Sporothrix_schenckii_S_SsMS1_GCA_002941045.out:masked repeats: 1,952,049 bp (5.98%)
Sporothrix_thermara_S_CBS139747_GCA_XXXXXXXXX.out:masked repeats: 2,829,844 bp (8.61%)
Sporothrix_thermara_S_bin10_GCA_XXXXXXXXX.out:masked repeats: 1,977,546 bp (7.22%)
Sporothrix_variecibatus_S_CBS121960_GCA_016097105.out:masked repeats: 2,682,515 bp (6.88%)
Sporothrix_variecibatus_S_CBS121961_GCA_021396255.out:masked repeats: 2,659,277 bp (6.84%)

