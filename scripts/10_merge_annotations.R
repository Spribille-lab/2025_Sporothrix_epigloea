#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

input_file <- args[1]
input_file <- "Sporothrix_thermara_S_CBS139747MAG_GCA_XXXXXXXXX"
input_file_short <- stringr::str_replace(input_file, "(\\S+)(_S)(_\\S+)(_GCA\\S+)", "\\1\\3")

# in the command line run as for example:
# Rscript get_secreted_proteases_from_deeploc.R ~/Documents/2022_04_eukaryote_annotation/08_secretome/ Sporothrix_epigloea_S_TF4163_GCA_94403644

## This assumes a fixed dir structure
# If you want you can create your dir structure
## dir.create("dir1", showWarnings = FALSE)
## dir.create("dir1/subdir1_1", showWarnings = FALSE)

library(seqRFLP)
library(tidyverse)
library(stringr)

setwd("~/Documents/2023_02_Sporothrix/10_annotations/")

annotations <- read.delim(paste("annotation_tables_funannotate/",input_file_short,".annotations.txt", sep=""), header = T, na.strings='')
colnames(annotations) = c("Gene_ID",
                         "Transcript_ID",
                         "Feature",
                         "Contig",
                         "Start",
                         "Stop",
                         "Strand",
                         "Name",
                         "Product",
                         "Alias_Synonyms",
                         "EC_number",
                         "BUSCO",
                         "PFAM",
                         "InterPro",
                         "EggNog",
                         "COG",
                         "GO_Terms",
                         "Secreted",
                         "Membrane",
                         "Protease",
                         "CAZyme",
                         "antiSMASH",
                         "Notes",
                         "gDNA",
                         "mRNA",
                         "CDS_transcript",
                         "Translation")

# read in the output from DeepLoc 2.0
deeploc <- read.csv(paste("DeepLoc2/",input_file,".deeploc2.csv", sep="")) %>%
  dplyr::mutate(Gene_ID = str_replace(Protein_ID, "(\\S*)(-T1)", "\\1")) %>%
  dplyr::relocate(Gene_ID) %>%
  dplyr::select(Gene_ID, Localizations, Signals)
colnames(deeploc) = c("Gene_ID", "Localizations_DeepLoc2", "Signals_DeepLoc2")

# merge DeepLoc and funannotate annotations

annotations2 <- dplyr::left_join(annotations, deeploc) %>%
  select(Gene_ID,
         Transcript_ID,
         Feature,
         Contig,
         Start,
         Stop,
         Strand,
         Name,
         Product,
         Alias_Synonyms,
         EC_number,
         BUSCO,
         PFAM,
         InterPro,
         EggNog,
         COG,
         GO_Terms,
         Secreted,
         Membrane,
         Localizations_DeepLoc2,
         Signals_DeepLoc2,
         Protease,
         CAZyme,
         antiSMASH,
         Notes,
         gDNA,
         mRNA,
         CDS_transcript,
         Translation)


# read in the output from OrthoFinder

orthogroups <- read.delim(paste("../09_orthology/OrthoFinder/Results_Aug07/Custom_results/",input_file_short,".orthogroups.txt", sep=""), header = T, na.strings='') %>%
  dplyr::mutate(Gene_ID = str_replace(Protein_ID, "(\\S*)(-T1)", "\\1")) %>%
  dplyr::select(Gene_ID, Orthogroup)

annotations3 <- dplyr::left_join(annotations2, orthogroups) %>%
  select(Gene_ID,
         Transcript_ID,
         Feature,
         Contig,
         Start,
         Stop,
         Strand,
         Name,
         Orthogroup,
         Product,
         Alias_Synonyms,
         EC_number,
         BUSCO,
         PFAM,
         InterPro,
         EggNog,
         COG,
         GO_Terms,
         Secreted,
         Membrane,
         Localizations_DeepLoc2,
         Signals_DeepLoc2,
         Protease,
         CAZyme,
         antiSMASH,
         Notes,
         gDNA,
         mRNA,
         CDS_transcript,
         Translation)

write.table(annotations3, paste("annotation_tables_merged/",input_file,".annotations.merged.txt", sep=""), sep='\t', quote = F, row.names = F)


##################################
#filter gene models labelled as a secreted protease

#protease <- annotations3[!is.na(annotations3$Protease),c('Gene_ID', 'Localizations_DeepLoc2','Translation')] %>%
#  dplyr::filter(Localizations_DeepLoc2 == "Extracellular") %>%
#  select(Gene_ID, Translation)

#write them as fasta
#protease_faa<-dataframe2fas(protease,file = paste("fastas/",input_file,".secretedprotease.faa", sep=""))


# filter rows with the proteinID found as secreted by DeepLoc
#secreted <- annotations3[annotations3$Gene_ID %in% deeploc$Gene_ID, c("Gene_ID","Translation")]
#deeploc_faa <- dataframe2fas(secreted,file = paste("secreted/",input_file,"_secreted.faa", sep=""))

# filter proteases that are secreted
#genome_secreted_proteases <- secreted[secreted$Gene_ID %in% protease$Gene_ID, c("Gene_ID","Translation")]
#genome_secreted_proteases_faa <- dataframe2fas(genome_secreted_proteases,file = paste("secreted_proteases_fastas/",input_file,"_secreted_proteases.faa", sep=""))

#protease_dataset_secreted <- protease_full[protease_full$Gene_ID %in% secreted$Gene_ID,]
#write.table(protease_dataset_secreted, paste("secreted_proteases_tables/",input_file,"_secreted_proteases.txt", sep=""), sep='\t', quote = F, row.names = F)

