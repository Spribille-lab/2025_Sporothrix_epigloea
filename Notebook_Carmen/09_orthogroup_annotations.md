
# Ran through dbCAN3 online server
2023.09.27
CAZyme domain boundaries of extracted orthogroups were annotated with the dbCAN3 webserver (Zheng et al., 2023) using a HMMER search of the dbCAN HMMdb release 12.0 containing 783 CAZyme HMMs.  Substrate annotations of CAZymes were enabled with a HMMER search of the CAZyme subfamily HMM database dbCAN-sub (released 2022.08.25).  For both, the HMMER thresholds were set to E-value < 1*10-15 and coverage > 0.35. Lastly, a DIAMOND alignment-based search of the Carbohydrate-Active enZYmes (CAZy) Database (CAZyDB; released on 2023.07.26) was performed with an E-value threshold < 1*10-102.  CAZyme proteins identified by two or three of the above tools were retained in further analyses.  


# annotate orthogroups gained and lost based on KEGG Orthology and hidden Markov models
- using KofamScan version 1.3.0
- 2023.09.27
- deBary

```
cd /data/ccallen/2023_02_Sporothrix/09_orthology/Results_Aug07/gained_and_lost
/data/ccallen/bin/kofam_scan/exec_annotation -o KEGG/gained.ko.txt epigloea_orthogroups_gained.faa -f detail-tsv &>> kofam_scan.gained.out
/data/ccallen/bin/kofam_scan/exec_annotation -o KEGG/gained.kegg.mapper.txt epigloea_orthogroups_gained.faa -f mapper &>> kofam_scan.gained.out

/data/ccallen/bin/kofam_scan/exec_annotation -o KEGG/lost.ko.txt epigloea_orthogroups_lost.faa -f detail-tsv &>> kofam_scan.lost.out
/data/ccallen/bin/kofam_scan/exec_annotation -o KEGG/lost.kegg.mapper.txt epigloea_orthogroups_lost.faa -f mapper &>> kofam_scan.lost.out
```

* move to local machine

