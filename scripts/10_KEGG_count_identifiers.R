#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

genome <- args[1]

## This assumes a fixed dir structure

library(tidyverse)
library(dplyr)

#genome = "Sporothrix_bragantina_S_CBS47491_GCA_XXXXXXXXX"

setwd(paste0("~/Documents/2023_02_Sporothrix/10_annotations/kegg/kegg_output/",genome, sep=""))
getwd()
kegg_mapper <- read.delim(paste(genome,".kegg.mapper.txt", sep=""), header = F, na.strings='') %>%
  select(-V1) %>%
  filter(!is.na(V2)) %>%
  arrange(V2)

kegg_mapper_count <- kegg_mapper %>% count(V2)
colnames(kegg_mapper_count) <- c(paste(genome), paste("n_",genome, sep=""))
write.table(kegg_mapper_count, paste(genome,".count_KEGG.txt", sep=""), sep='\t', quote = F, row.names = F)
