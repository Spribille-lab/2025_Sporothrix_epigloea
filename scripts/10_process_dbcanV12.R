

## 1. misc

library(tidyverse)
library(phylotools)


setwd("~/Documents/2023_02_Sporothrix/10_annotations/dbcan4/")

## 2. read data

annotated_genomes<-read.delim("../genomes_dbcan.txt", header=T)

## 4. make a function to process dbcan output

#target_genome = "Sporothrix_bragantina_S_CBS47491_GCA_XXXXXXXXX" #to test function

#overview_filename <- paste0(target_genome,"_dbcan/overview.txt")
#overview <-read.delim(overview_filename) %>% filter(X.ofTools > 1)

# hits from hmmer and diamond ##### I didn't do this in the end #######
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

#this one only looks at hmmer results and doesn't give the regions########
process_dbcan_hmmer<-function(target_genome){
  
  hmmer_filename<-paste0(target_genome,"_dbcan/hmmer.out")
  hmmer<-read.delim(hmmer_filename)
  
  hmmer_filtered <- hmmer %>% filter(E.Value<1e-20,Coverage>0.3)
  hmmer_filtered$family <- str_replace(hmmer_filtered$HMM.Profile,".hmm","")
  
  #Filter for AA14s
  combined <- hmmer_filtered %>% select(Gene.ID, family) #%>% filter(family == "AA14")
  combined$genome<-target_genome

  return(combined)
}

# for fungi, use E-value < 1e-17 and coverage > 0.45.

# hmmer hits with regions.  This is the one that I used. ###########
process_dbcan_hmmer_regions<-function(target_genome){
  
  hmmer_filename<-paste0(target_genome,"_dbcan/hmmer.out")
  hmmer<-read.delim(hmmer_filename)
  
  hmmer_filtered <- hmmer %>% filter(E.Value<1e-17,Coverage>0.45)
  hmmer_filtered$family <- str_replace(hmmer_filtered$HMM.Profile,".hmm","")
  
  combined <- hmmer_filtered %>% select(Gene.ID, family, Gene.Start, Gene.End)
  combined$genome<-target_genome
  
  return(combined)
}


## 5. apply the function to all genomes
l<-lapply(annotated_genomes$genome,process_dbcan_hmmer_regions)
cazy_combined <- do.call(rbind,l)
cazy_combined$class <- substr(cazy_combined$family,1,2)

cazy_combined <- cazy_combined %>%
  #filter(family == "AA9" | family == "AA11" | family == "AA13" | family == "AA14" | family == "AA15" | family == "AA16" | family == "AA17") %>%
  unite("region", c("Gene.ID", "Gene.Start", "Gene.End"), sep = "_", remove = FALSE) %>%
  mutate(region = str_replace(region, "(\\S+-T1)_(\\d+)_(\\d+)","\\1:\\2-\\3")) %>%
  unite("tip_labels", c("Gene.ID", "Gene.Start", "Gene.End", "family"), sep = "_", remove = FALSE) %>%
  mutate(tip_labels = str_replace(tip_labels, "(\\S+-T1)_(\\d+)_(\\d+)(_AA\\d+)","\\1:\\2-\\3\\4")) %>%
  left_join(annotated_genomes, by="genome")

lpmo_combined <- cazy_combined %>%
  filter(family == "AA9" | family == "AA11" | family == "AA13" | family == "AA14" | family == "AA15" | family == "AA16" | family == "AA17")

lpmo_combined_AA14 <- lpmo_combined %>%
  filter(family == "AA14")

lpmo_combined_AA14_sordariomycetes <- lpmo_combined %>%
  filter(family == "AA14") %>%
  filter(tax_class == "Sordariomycetes")

# write a table with the gene name, start and end positions, region (how each tip is named in iqtree), tip_labels (how I want to name the tips)
write.table(lpmo_combined, "summarized_outputs/LPMO_gene_assignments.txt", sep='\t',quote = F, row.names = F, col.names = T)
write.table(lpmo_combined_AA14, "summarized_outputs/AA14_gene_assignments.txt", sep='\t',quote = F, row.names = F, col.names = T)
write.table(lpmo_combined_AA14_sordariomycetes, "summarized_outputs/AA14_Sordariomycetes_gene_assignments.txt", sep='\t',quote = F, row.names = F, col.names = T)

lpmo_combined %>% distinct(genome)

# make a table to extract sequences from fasta using the gene name and the start and end positions
extract_sequence_list <- lpmo_combined %>%
  mutate(fasta = str_replace(genome, "(.+)","\\1.proteins.fa")) %>%
  select(genome, Gene.ID, Gene.Start, Gene.End)
write.table(extract_sequence_list, "summarized_outputs/lpmo_coordinates.txt", sep='\t',quote = F, row.names = F, col.names = F)
  
extract_sequence_list_AA14 <- lpmo_combined_AA14 %>%
  mutate(fasta = str_replace(genome, "(.+)","\\1.proteins.fa")) %>%
  select(genome, Gene.ID, Gene.Start, Gene.End)
write.table(extract_sequence_list_AA14, "summarized_outputs/AA14_coordinates.txt", sep='\t',quote = F, row.names = F, col.names = F)

extract_sequence_list_AA14_Sordariomycetes <- lpmo_combined_AA14_sordariomycetes %>%
  mutate(fasta = str_replace(genome, "(.+)","\\1.proteins.fa")) %>%
  select(genome, Gene.ID, Gene.Start, Gene.End)
write.table(extract_sequence_list_AA14_Sordariomycetes, "summarized_outputs/AA14_Sordariomycetes_coordinates.txt", sep='\t',quote = F, row.names = F, col.names = F)


# make one fasta file with AA14/LPMO sequences of all genomes
#```
#cd /data/ccallen/2023_02_Sporothrix/10_annotations/dbcan4
#while IFS=$'\t' read -r v1 v2 v3 v4; do samtools faidx ../protein_fasta_dbcan/"$v1".proteins.fa $v2:$v3-$v4 >> summarized_outputs/allgenomes.cazy_aa14.fa; done < summarized_outputs/AA14_coordinates.txt
#while IFS=$'\t' read -r v1 v2 v3 v4; do samtools faidx ../protein_fasta_dbcan/"$v1".proteins.fa $v2:$v3-$v4 >> summarized_outputs/Sordariomycetes.cazy_aa14.fa; done < summarized_outputs/AA14_Sordariomycetes_coordinates.txt
#while IFS=$'\t' read -r v1 v2 v3 v4; do samtools faidx ../protein_fasta_dbcan/"$v1".proteins.fa $v2:$v3-$v4 >> summarized_outputs/allgenomes.cazy_lpmo.fa; done < summarized_outputs/lpmo_coordinates.txt
#```

ref_table <- lpmo_combined %>%
  select(region, tip_labels)
phylotools::rename.fasta(infile = "summarized_outputs/allgenomes.cazy_aa14.fa", ref_table, outfile = "summarized_outputs/allgenomes.cazy_aa14.renamed.fa")
phylotools::rename.fasta(infile = "summarized_outputs/allgenomes.cazy_lpmo.fa", ref_table, outfile = "summarized_outputs/allgenomes.cazy_lpmo.renamed.fa")



aa14_fasta <- read.fasta("summarized_outputs/allgenomes.cazy_aa14.renamed.fa")




## 6. process the result
cazy_summarized <- cazy_combined %>% group_by(genome,family) %>% summarize(n=n()) 


#human-readable table
hr_table <- cazy_summarized %>% pivot_wider(names_from=family,values_from=n,values_fill = 0) %>% 
  left_join(annotated_genomes) %>%
  relocate(tax_class, .after = genome)

write.table(hr_table, "summarized_outputs/cazy_summarized.txt", sep='\t',quote = F, row.names = F, col.names = T)

## 7. make a list of genes for each genome -------

# make a function to pull names of AA9, 11, 13 and 14 genes for each genome

pull_AAs <-function(target_genome){
  AAgenes <- cazy_combined %>% filter(genome == target_genome) %>% pull(Gene.ID)
  return(AAgenes)
}

# apply the function to all genomes

l <- sapply(annotated_genomes$genome, pull_AAs)







##### make a heatmap of CAZymes ######


library(data.table)
library(ComplexHeatmap)


###################tree#############
library("ape")
library('dendextend')
library("DECIPHER")


ML_tree <- ape::read.tree("~/Documents/2023_02_Sporothrix/09_orthology/OrthoFinder/Results_Aug03/Species_Tree/iqtree2/sporothrix_placement.rerooted.contree")
ML_tree$tip.label
ML_tree <- drop.tip(ML_tree, 25)
str(ML_tree)

ML_dend <- ReadDendrogram(textConnection(write.tree(ML_tree)))
str(ML_dend)
length(ML_tree$tip.label)
ML_dend %>% nleaves

order.dendrogram(ML_dend)
#order.dendrogram(ML_dend) <- 1:length(ML_tree$tip.label)
order.dendrogram(ML_dend) <- 1:nleaves(ML_dend)
order.dendrogram(ML_dend)

#make a dataframe of genomes.  I use this dataframe as an anchor to keep the dendrogram and the data aligned.
genome <- unlist(ML_tree$tip.label) 
genome_df <- data.frame(genome)


# make a list of the CAZyme families and their putative functions

CAZyme_activity <- read.delim("~/Documents/2022_04_eukaryote_annotation/06_dbcan/CAZyme_activity.txt", header = T, na.strings='')
CAZyme_activity <- CAZyme_activity %>% select(General_Activity, CAZyme_families, origin)
colnames(CAZyme_activity) = c("General_Activity", "CAZyme_family", "origin")



# make matrix for full heatmap
cazyme_matrix <- hr_table %>% as.data.frame() %>% select(-tax_class) %>%
  mutate(genome = str_replace(genome, "(\\S+)(_\\S+)_S(_\\S+)(_GCA_\\S+)","\\1\\2\\3"))
cazyme_matrix$genome[cazyme_matrix$genome == 'Sporothrix_thermara_CBS139747'] <- 'Sporothrix_thermara_CBS139747MAG'
str(cazyme_matrix)
cazyme_matrix <- left_join(genome_df, cazyme_matrix)
row.names(cazyme_matrix) = cazyme_matrix$genome
cazyme_matrix <- cazyme_matrix %>% dplyr::select(-genome)

cazyme_matrix_t <- as.data.frame(t(cazyme_matrix))
cazyme_matrix_t <- cazyme_matrix_t %>% rownames_to_column(var = "CAZyme_family")

cazyme_matrix_t <- cazyme_matrix_t %>%
  dplyr::left_join(CAZyme_activity, by = "CAZyme_family") %>%
  dplyr::arrange(General_Activity) %>%
  dplyr::relocate(General_Activity, .after = CAZyme_family) %>%
  dplyr::relocate(origin, .after = CAZyme_family)

#######make a table for CAFE########
hr_table_for_CAFE <- cazyme_matrix_t %>% select (-origin) %>% relocate(General_Activity)
colnames(hr_table_for_CAFE)[1] <- "Desc"
colnames(hr_table_for_CAFE)[2] <- "Family ID"

write.table(hr_table_for_CAFE, "summarized_outputs/cazy_summarized_CAFE.txt", sep='\t',quote = F, row.names = F, col.names = T)


# back to heatmap
row.names(cazyme_matrix_t) = cazyme_matrix_t$CAZyme_family
cazyme_matrix_t <- cazyme_matrix_t %>% select(-CAZyme_family)



cazyme_matrix <- as.matrix(cazyme_matrix)
#cazyme_matrix <- apply(cazyme_matrix, 8236, as.numeric)
class(cazyme_matrix)

CAZyme_activity_matrix <- cazyme_matrix_t %>%
  select(origin, General_Activity) %>%
  as.matrix()

cazyme_matrix_t <- cazyme_matrix_t %>%
  select(-origin, -General_Activity) %>%
  as.matrix()



############### vertical heatmap #############

library(circlize)
col_fun = colorRamp2(c(0, max(unlist(cazyme_matrix_t))), c("white", "darkorchid"))
max(unlist(cazyme_matrix))

family_annotation <- rowAnnotation(
  "General_activity" = anno_text(CAZyme_activity_matrix[,2], gp = gpar(fontsize = 6)),
  "origin" = anno_text(CAZyme_activity_matrix[,1], gp = gpar(fontsize = 6)))

cazyme_heatmap <- Heatmap(cazyme_matrix_t,
                          name = "n",
                          row_order = order(CAZyme_activity_matrix[, 1]),
                          column_names_gp = grid::gpar(fontsize = 7),
                          row_names_gp = grid::gpar(fontsize = 6),
                          row_names_side = "left",
                          #cluster_rows = FALSE,
                          col = col_fun,
                          width = unit(9, "cm"),
                          height = unit(70, "cm"),
                          right_annotation = family_annotation,
                          show_heatmap_legend = FALSE,
                          cell_fun = function(j, i, x, y, width, height, fill) {
                            grid.text(cazyme_matrix_t[i, j], x, y, gp = gpar(fontsize = 6))}
)
#draw(cazyme_heatmap)

pdf(file="../figures/sporothrix_allcazyme_heatmap_complexheatmap.pdf", width=10, height=30)

draw(cazyme_heatmap)
dev.off()
