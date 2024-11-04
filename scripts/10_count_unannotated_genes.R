library(tidyverse)

setwd("~/Documents/2023_02_Sporothrix/10_annotations/")


file_list <- list.files(path = "annotation_tables_merged/")
target_file = "Leptographium_lundbergii_S_CBS138716_GCA_001455505.annotations.merged.txt"

count_unannotated_genes <- function(target_file) {
  
  filename <- paste0("annotation_tables_merged/",target_file)
  file <- read.delim(filename)
  
  genes <- nrow(file)
  
  file_NA <- file %>%
    select(Gene_ID, KEGG_ID, PFAM, InterPro, Protease, CAZyme) %>%
    filter(is.na(KEGG_ID)) %>%
    filter(is.na(PFAM)) %>%
    filter(is.na(InterPro)) %>%
    filter(is.na(Protease)) %>%
    filter(is.na(CAZyme))
  
  unannotated <- nrow(file_NA)
  unannotated <- as.data.frame(unannotated)
  
  unannotated$genes <- genes
  
  unannotated$genome <- target_file
  unannotated <- unannotated %>%
    mutate(genome = str_replace(genome, "(\\S+)(.annotations.merged.txt)","\\1"))
  
  return(unannotated)
}

l<-lapply(file_list,count_unannotated_genes)
unannotated_genes <- do.call(rbind,l)
write.table(unannotated_genes, "unannotated_genes.txt", sep='\t',quote = F, row.names = F, col.names = T)
