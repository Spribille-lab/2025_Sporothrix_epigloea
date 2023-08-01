

## 1. misc

library(tidyverse)

setwd("~/Documents/2023_02_Sporothrix/10_annotations/dbcan3/")

## 2. read data

annotated_genomes<-read.delim("../genomes.txt", header=T)

## 4. make a function to process dbcan output

target_genome = "Sporothrix_bragantina_S_CBS47491_GCA_XXXXXXXXX" #to test function

overview_filename <- paste0(target_genome,"_dbcan/overview.txt")
overview <-read.delim(overview_filename) %>% filter(X.ofTools > 1)

# hits from hmmer and diamond
process_dbcan_2<-function(target_genome){
  
  hmmer_filename<-paste0(target_genome,"_dbcan/hmmer.out")
  hmmer<-read.delim(hmmer_filename)
  
  hmmer_filtered <- hmmer %>% filter(E.Value<1e-20,Coverage>0.3)
  hmmer_filtered$subfamily <- str_replace(hmmer_filtered$HMM.Profile,".hmm","")
  
  diamond_filename <- paste0(target_genome,"_dbcan/diamond.out")
  diamond <- read.delim(diamond_filename)
  
  diamond_filtered <- diamond %>% filter(E.Value<1e-20,X..Identical>30)
  diamond_filtered <- diamond_filtered %>% mutate(subfamily = str_replace(CAZy.ID, "^(.+?)\\|","")) %>% 
    separate_rows(subfamily,sep = "\\|", convert = FALSE)
  
  dbsub_filename <- paste0(target_genome,"_dbcan/dbsub.out")
  dbsub <- read.delim(dbsub_filename)
  
  dbsub_filtered <- dbsub %>% filter(E.Value<1e-20,Coverage>0.3)
  dbsub_filtered <- dbsub_filtered %>% mutate(family = str_replace(dbCAN.subfam, "\\_\\S*","")) %>%
    select(Gene.ID, family, Substrate)
  
  #combine diamond and hmmer results and keep only those with matches.  Filter for AA14s
  combined <- diamond_filtered %>% select(Gene.ID, subfamily) %>% inner_join(hmmer_filtered) %>%
    select(Gene.ID, subfamily) # %>% filter(family == "AA14")
  combined$genome <- target_genome
  return(combined)
}

#this one only looks at hmmer results
process_dbcan_hmmer<-function(target_genome){
  
  hmmer_filename<-paste0("dbcan_output/",target_genome,"_dbcan/hmmer.out")
  hmmer<-read.delim(hmmer_filename)
  
  hmmer_filtered <- hmmer %>% filter(E.Value<1e-20,Coverage>0.3)
  hmmer_filtered$family <- str_replace(hmmer_filtered$HMM.Profile,".hmm","")
  
  #Filter for AA14s
  combined <- hmmer_filtered %>% select(Gene.ID, family) #%>% filter(family == "AA14")
  combined$genome<-target_genome

  return(combined)
}

# hmmer hits with regions
process_dbcan_hmmer_regions<-function(target_genome){
  
  hmmer_filename<-paste0("dbcan_output/",target_genome,"_dbcan/hmmer.out")
  hmmer<-read.delim(hmmer_filename)
  
  hmmer_filtered <- hmmer %>% filter(E.Value<1e-20,Coverage>0.3)
  hmmer_filtered$family <- str_replace(hmmer_filtered$HMM.Profile,".hmm","")
  
  #Filter for AA14s
  combined <- hmmer_filtered %>% select(Gene.ID, family, Gene.Start, Gene.End) #%>% filter(family == "AA14")
  combined$genome<-target_genome
  
  return(combined)
}


## 5. apply the function to all genomes
l<-lapply(annotated_genomes$genome,process_dbcan_hmmer_regions)
cazy_combined<-do.call(rbind,l)
cazy_combined$class<-substr(cazy_combined$family,1,2)
cazy_combined <- cazy_combined %>%
  filter(family == "AA9" | family == "AA11" | family == "AA13" | family == "AA14") %>%
  unite("region", c("Gene.ID", "Gene.Start", "Gene.End"), sep = "_", remove = FALSE) %>%
  mutate(region = str_replace(region, "(\\S+-T1)_(\\d+)_(\\d+)","\\1_\\2-\\3")) %>%
  unite("tip_labels", c("Gene.ID", "Gene.Start", "Gene.End", "family"), sep = "_", remove = FALSE) %>%
  mutate(tip_labels = str_replace(tip_labels, "(\\S+-T1)_(\\d+)_(\\d+)(_AA\\d+)","\\1:\\2-\\3\\4"))


# write a table with the gene name, start and end positions, region (how each tip is named in iqtree), tip_labels (how I want to name the tips)
write.table(cazy_combined, "summarized_outputs/AA14_gene_assignments.txt", sep='\t',quote = F, row.names = F, col.names = T)

cazy_combined %>% distinct(genome)

# make a table to extract sequences from fasta using the gene name and the start and end positions
extract_sequence_list <- cazy_combined %>%
  mutate(fasta = str_replace(genome, "(.+)","\\1.proteins.fa")) %>%
  select(genome, Gene.ID, Gene.Start, Gene.End)
write.table(extract_sequence_list, "summarized_outputs/lpmo_coordinates.txt", sep='\t',quote = F, row.names = F, col.names = F)
  




## 6. process the result
cazy_summarized<-cazy_combined %>% group_by(genome,family) %>% summarize(n=n()) 


#human-readable table
hr_table<-cazy_summarized %>% pivot_wider(names_from=family,values_from=n,values_fill = 0) %>% 
  left_join(annotated_genomes) 

## 7. make a list of genes for each genome -------

# make a function to pull names of AA9, 11, 13 and 14 genes for each genome

pull_AAs <-function(target_genome){
  AAgenes <- cazy_combined %>% filter(genome == target_genome) %>% pull(Gene.ID)
  return(AAgenes)
}

# apply the function to all genomes

l <- sapply(annotated_genomes$genome, pull_AAs)





