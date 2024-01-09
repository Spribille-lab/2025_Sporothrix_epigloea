# check completeness of each CBS139747 bin using busco

* 2023.06.15
* in debary


# copy all concoct bins into a single directory and renamed

```
cd /data/ccallen/2023_02_Sporothrix/06_binning
for sample in $(cat ../samples.txt)
do
mkdir all_concoct_bins/"$sample"
cd "$sample"/concoct_output/fasta_bins
for file in *.fa
do
echo $file
cp $file ../../../all_concoct_bins/"$sample"
mv ../../../all_concoct_bins/"$sample"/$file ../../../all_concoct_bins/"$sample"_$file
done
cd ../../..
rm -r all_concoct_bins/"$sample"
done
```




Run busco on the CBS139747 bins

```
cd /data/ccallen/2023_02_Sporothrix/07_bin_QC
for path in $(ls bins/*.fa)
do
file=${path##*/}
bin=${file%.*}
echo $bin  
docker run -u $(id -u) -v $(pwd):/busco_wd ezlabgva/busco:v5.4.7_cv1 busco -i "$path" -o "$bin" -m genome --cpu 8 -f
done
```

compile report
```
cd /data/ccallen/2023_02_Sporothrix/07_bin_QC
for filename in */*.txt
do
echo -n  $filename >> busco_report.txt;
sed -n -e '/^\tC:/p' $filename >> busco_report.txt
done
```

CBS139747_0/short_summary.specific.eukaryota_odb10.CBS139747_0.txt	C:2.0%[S:2.0%,D:0.0%],F:0.8%,M:97.2%,n:255	   
CBS139747_10/short_summary.generic.eukaryota_odb10.CBS139747_10.txt	C:92.2%[S:92.2%,D:0.0%],F:3.9%,M:3.9%,n:255	   
CBS139747_10/short_summary.specific.sordariomycetes_odb10.CBS139747_10.txt	C:87.7%[S:87.5%,D:0.2%],F:3.5%,M:8.8%,n:3817	   
CBS139747_2/short_summary.specific.eukaryota_odb10.CBS139747_2.txt	C:0.8%[S:0.8%,D:0.0%],F:0.0%,M:99.2%,n:255	   
CBS139747_23/short_summary.generic.eukaryota_odb10.CBS139747_23.txt	C:98.0%[S:98.0%,D:0.0%],F:0.8%,M:1.2%,n:255	   
CBS139747_23/short_summary.specific.sordariomycetes_odb10.CBS139747_23.txt	C:95.4%[S:95.2%,D:0.2%],F:1.2%,M:3.4%,n:3817	   
CBS139747_24/short_summary.specific.eukaryota_odb10.CBS139747_24.txt	C:2.0%[S:2.0%,D:0.0%],F:0.4%,M:97.6%,n:255	   
CBS139747_34/short_summary.generic.bacteria_odb10.CBS139747_34.txt	C:0.8%[S:0.8%,D:0.0%],F:0.0%,M:99.2%,n:124	   
CBS139747_34/short_summary.specific.mollicutes_odb10.CBS139747_34.txt	C:0.0%[S:0.0%,D:0.0%],F:0.7%,M:99.3%,n:151	   
CBS139747_8/short_summary.specific.archaea_odb10.CBS139747_8.txt	C:0.5%[S:0.5%,D:0.0%],F:0.0%,M:99.5%,n:194

```
docker run -u $(id -u) -v $(pwd):/busco_wd ezlabgva/busco:v5.4.7_cv1 busco -i bins/CBS139747_1_2_3_21_23_26_34.fa -o CBS139747_1_2_3_21_23_26_34 -m genome --cpu 8 -f
```

CBS139747_1_2_3_21_23_26_34/short_summary.generic.eukaryota_odb10.CBS139747_1_2_3_21_23_26_34.txt	C:98.8%[S:98.8%,D:0.0%],F:0.8%,M:0.4%,n:255	   
CBS139747_1_2_3_21_23_26_34/short_summary.specific.sordariomycetes_odb10.CBS139747_1_2_3_21_23_26_34.txt	C:95.6%[S:95.4%,D:0.2%],F:1.2%,M:3.2%,n:3817	   
CBS139747_1_2_3_9_21_23_26_33_34/short_summary.specific.sordariomycetes_odb10.CBS139747_1_2_3_9_21_23_26_33_34.txt	C:95.6%[S:95.4%,D:0.2%],F:1.2%,M:3.2%,n:3817	





* short_summary.specific.sordariomycetes_odb10.CBS139747_10.txt
		***** Results: *****

	C:87.7%[S:87.5%,D:0.2%],F:3.5%,M:8.8%,n:3817	   
	3344	Complete BUSCOs (C)			   
	3338	Complete and single-copy BUSCOs (S)	   
	6	Complete and duplicated BUSCOs (D)	   
	133	Fragmented BUSCOs (F)			   
	340	Missing BUSCOs (M)			   
	3817	Total BUSCO groups searched		   

Assembly Statistics:
	1048	Number of scaffolds
	1048	Number of contigs
	27371975	Total length
	0.000%	Percent gaps
	37 KB	Scaffold N50
	37 KB	Contigs N50

* short_summary.specific.sordariomycetes_odb10.CBS139747_23.txt
	***** Results: *****

	C:95.4%[S:95.2%,D:0.2%],F:1.2%,M:3.4%,n:3817	   
	3639	Complete BUSCOs (C)			   
	3633	Complete and single-copy BUSCOs (S)	   
	6	Complete and duplicated BUSCOs (D)	   
	47	Fragmented BUSCOs (F)			   
	131	Missing BUSCOs (M)			   
	3817	Total BUSCO groups searched		   

Assembly Statistics:
	192	Number of scaffolds
	192	Number of contigs
	32530678	Total length
	0.000%	Percent gaps
	276 KB	Scaffold N50
	276 KB	Contigs N50


* short_summary.specific.sordariomycetes_odb10.CBS139747_1_2_3_21_23_26_34.txt #this is the one that we used!!
		***** Results: *****

	C:95.6%[S:95.4%,D:0.2%],F:1.2%,M:3.2%,n:3817	   
	3646	Complete BUSCOs (C)			   
	3640	Complete and single-copy BUSCOs (S)	   
	6	Complete and duplicated BUSCOs (D)	   
	47	Fragmented BUSCOs (F)			   
	124	Missing BUSCOs (M)			   
	3817	Total BUSCO groups searched		   

Assembly Statistics:
	205	Number of scaffolds
	205	Number of contigs
	32873637	Total length
	0.000%	Percent gaps
	276 KB	Scaffold N50
	276 KB	Contigs N50

* short_summary.specific.sordariomycetes_odb10.CBS139747_1_2_3_9_21_23_26_33_34.txt

		***** Results: *****

	C:95.6%[S:95.4%,D:0.2%],F:1.2%,M:3.2%,n:3817	   
	3646	Complete BUSCOs (C)			   
	3640	Complete and single-copy BUSCOs (S)	   
	6	Complete and duplicated BUSCOs (D)	   
	47	Fragmented BUSCOs (F)			   
	124	Missing BUSCOs (M)			   
	3817	Total BUSCO groups searched		   

Assembly Statistics:
	207	Number of scaffolds
	207	Number of contigs
	32885059	Total length
	0.000%	Percent gaps
	276 KB	Scaffold N50
	276 KB	Contigs N50


