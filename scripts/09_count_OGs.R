library(tidyverse)
library(phylotools)
setwd("~/Documents/2023_02_Sporothrix/09_orthology/OrthoFinder/Results_Aug07/gained_and_lost/")
lost_an <- read.delim("epigloea_orthogroups_lost_annotations.txt", stringsAsFactors = F)
lost_cazymes_an <- lost_an[!is.na(lost_an$CAZyme),]
lost_cazymes_og <- lost_cazymes_an %>% distinct(Orthogroup)
lost_cazymes_fam <- lost_cazymes_an %>% distinct(CAZyme) %>% separate_rows(CAZyme,sep = ";", convert = FALSE) %>% distinct(CAZyme)


GH <- filter(lost_cazymes_fam, grepl("GH",CAZyme))
GH$CAZyme


lost_og <- lost_an %>% distinct(Orthogroup)

##proteases
lost_proteases_an <- lost_an[!is.na(lost_an$Protease),]
lost_proteases_og <- lost_proteases_an %>% distinct(Orthogroup)
lost_proteases_fam <- lost_proteases_an %>% distinct(Protease) %>% sort(Protease)

sort(lost_proteases_fam$Protease)

S <- filter(lost_proteases_fam, grepl("S",Protease))
S$Protease

#######transporters#######

MFS <- filter(lost_an, grepl("Major facilitator superfamily|MFS", InterPro))
MFS_og <- MFS %>% distinct(Orthogroup)

PF00083 <- filter(lost_an, grepl("PF00083", PFAM))

#sugar_word <- filter(lost_an, grepl("Sugar transporter|sugar transporter", InterPro))
sugar_IPR <- filter(lost_an, grepl("IPR005829|IPR003663|IPR005828", InterPro))
# sugar_join <- sugar_word |> anti_join(sugar_IPR, join_by(Transcript_ID)) # test which are missing
#sugar_join <- PF00083 |> anti_join(sugar_IPR, join_by(Transcript_ID)) # test which are missing

sugar_transporters <- bind_rows(sugar_IPR, PF00083)
sugar_transporters_og <- sugar_transporters %>% distinct(Orthogroup)

non_sugar_MFS <- MFS %>% anti_join(sugar_transporters, join_by(Transcript_ID))
non_sugar_MFS_og <- non_sugar_MFS %>% distinct(Orthogroup)

## other stuff #####

others_temp <- lost_an[is.na(lost_an$CAZyme),]
others2 <- others_temp[is.na(others_temp$Protease),]

others3 <- filter(others2, !grepl("Major facilitator superfamily|MFS", InterPro))
others <- filter(others3, !grepl("PF00083", PFAM))
others_og  <- others %>% distinct(Orthogroup)

finger <- filter(others, grepl("zinc finger|Zn|SANT|DNA-binding|Transcription", InterPro))
finger_og  <- finger %>% distinct(Orthogroup)

Zn <- filter(others, grepl("Zn|zinc finger", InterPro))
Zn_og  <- Zn %>% distinct(Orthogroup)

SANT <- filter(others, grepl("SANT", InterPro))
SANT_og  <- SANT %>% distinct(Orthogroup)

BTB_POZ <- filter(others, grepl("IPR000210", InterPro))
BTB_POZ_og <- BTB_POZ %>% distinct(Orthogroup)

Ubiquitin <- filter(others, grepl("Ubiquitin|ubiquitin", InterPro)) 

signal_transduction <- filter(others, grepl("Signal transduction", COG)) 
signal_transduction_og <- signal_transduction %>% distinct(Orthogroup)

toxin  <- filter(others, grepl("toxin", InterPro)) 
toxin_og  <- toxin %>% distinct(Orthogroup)

nt_fasta_lost <- lost_an %>% select(Transcript_ID, gDNA)
colnames(nt_fasta_lost) <- c("seq.name","seq.text")
dat2fasta(nt_fasta_lost, outfile = "epigloea_orthogroups_lost.fna") # for antismash 7

gained_an <- read.delim("epigloea_orthogroups_gained_annotations.txt", stringsAsFactors = F)
nt_fasta_gained <- gained_an %>% select(Transcript_ID, gDNA)
colnames(nt_fasta_gained) <- c("seq.name","seq.text")
dat2fasta(nt_fasta_gained, outfile = "epigloea_orthogroups_gained.fna")
 

antismash <- read.delim("antismash_lost.txt", stringsAsFactors = F)
antismash_an <- left_join(antismash, lost_an, by="Transcript_ID")

