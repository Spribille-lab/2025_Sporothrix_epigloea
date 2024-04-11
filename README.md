# Gene-loss-in-Sporothrix-epigloea-shift-to-life-in-the-Tremella-fungal-heteropolymer-matrix

This repository contains scripts and intermediate results for the manuscript (Allen et al. 2024, <journal>)

## Abstract
The unique microhabitat of Sporothrix epigloea—the extracellular matrix of Tremella jelly fungi—is rich with acidic glucuronoxylomannan (GXM) polysaccharides.  GXM polysaccharides are medically relevant both as a pharmaceutical extracted from Tremella fruiting bodies and as a critical virulence factor involved in cryptococcosis caused by the tremelloid fungus Cryptococcus.  Using a comparative genomics approach, we surveyed the genome assemblies of Sporothrix species through the evolutionary lifestyle transition from soil or bark beetle associate to life on a GXM-rich fungal fruiting body.  In addition, we inferred protein orthogroups across core Sporothrix genomes to identify gene families that were either lost or derived in the divergence of S. epigloea.  We found that S. epigloea genomes were smaller than any other Sporothrix genome and reduced in carbohydrate active enzymes (CAZymes), proteases, biosynthetic gene clusters, and sugar transporters.  The suite of CAZyme families degrading both plant and fungal cell wall components were reduced in S. epigloea.  At the same time, a lytic polysaccharide monooxygenase (LPMO) with no previously established activity or substrate specificity, appears to have been uniquely acquired in S. epigloea.  This result calls for further investigation in the substrate activity of these enzymes and if they play a role in degrading GXM in Tremella fruiting bodies.


## Overview
```
project
├── README.md							# this doc; description of the repo and the project log
├── scripts							# all scripts generated for the analysis, with the exception of snakemake pipelines
├── 04_assemblies				
├── 08_GC_coverage						
├── 09_orthology
├── 10_annotations			
├── 11_LPMO_trees						
└── Notebook 						# temp log files for various parts of the analysis; the cleaned-up version of the same logs is in this file
└── results 							# results included in the manuscript
    ├── figures 						# figures
    └── tables 							# tables

```