library(tidyverse)
library(patchwork)
library(egg)
library(cowplot)
library("ComplexHeatmap")
library("ape")
library('dendextend')
library("DECIPHER")

sessionInfo(package = NULL)
#############################

# make a heatmap using complexheatmap that compares sporothrix species
setwd("~/Documents/2023_02_Sporothrix/10_annotations/")

ML_tree <- ape::read.tree("~/Documents/2023_02_Sporothrix/09_orthology/OrthoFinder/Restuls_Jul26/SpeciesTree_rooted_newick.txt")
ML_tree$tip.label

ML_dend <- ReadDendrogram(textConnection(write.tree(ML_tree)))
str(ML_dend)

#ML_dend <- DECIPHER::ReadDendrogram("~/Documents/2022_04_eukaryote_annotation/05_genomic_analysis_Sporothrix/ML_tree2.txt",keepRoot = T)
order.dendrogram(ML_dend)
order.dendrogram(ML_dend) <- 1:26 # change the last number to the number of tips.  I'm sure there's a fancy way to put the number of tips but I'm not that fancy.
order.dendrogram(ML_dend)

#make a dataframe of genomes.  I use this dataframe as an anchor to keep the dendrogram and the data aligned.
genome <- unlist(ML_tree$tip.label) 
genome_df <- data.frame(genome)

#make a dataframe of genome attributes
genome_attributes <- read.delim("~/Documents/2023_02_Sporothrix/10_annotations/genome_attributes.txt")
genome_attributes <- left_join(genome_df,genome_attributes)
row.names(genome_attributes) = genome_attributes$genome
genome_attributes_lifestyle <- genome_attributes %>%
  dplyr::select(-genome) %>%
  dplyr::select(mammal_pathogen, arthropod, wood, soil, infructescence, macrofungus)
str(genome_attributes_lifestyle)
genome_attributes_lifestyle_matrix <- as.matrix(genome_attributes_lifestyle)
str(genome_attributes_lifestyle_matrix)

#read in data
completeness <- read.delim("eukcc2/eukcc2.txt")
completeness <- left_join(genome_df,completeness)
completeness <- completeness %>% dplyr::select(-genome_assembly)
row.names(completeness) = completeness$genome
completeness <- completeness %>% dplyr::select(-genome)
completeness_matrix <- as.matrix(completeness)

genome_stats <-read.csv("funannotate/funannotate_compare_20230729/stats/genome.stats.summary.csv")
row.names(genome_stats) = genome_stats$X
genome_stats <- genome_stats %>% dplyr::select(-X)
genome_stats <- as.data.frame(t(genome_stats)) #transpose
rownames(genome_stats)
genome_stats <- genome_stats %>% add_column(genome = rownames(genome_stats), .before = "isolate") #make new column with the row names
genome_stats <- left_join(genome_df,genome_stats) 
row.names(genome_stats) = genome_stats$genome
genome_stats <- genome_stats %>% dplyr::select(-genome, -isolate, -locus_tag)
genome_stats <- transform(genome_stats,
                          Assembly_Size = as.numeric(Assembly_Size),
                          Average_Scaffold = as.numeric(Average_Scaffold),
                          Largest_Scaffold = as.numeric(Largest_Scaffold),
                          Num_Scaffolds = as.numeric(Num_Scaffolds),
                          Scaffold_N50  = as.numeric(Scaffold_N50),
                          Percent_GC = as.numeric(Percent_GC),
                          Num_Genes = as.numeric(Num_Genes),
                          Num_Proteins = as.numeric(Num_Proteins),
                          Num_tRNA = as.numeric(Num_tRNA),
                          Unique_Proteins = as.numeric(Unique_Proteins),
                          Prots_atleast_1_ortholog = as.numeric(Prots_atleast_1_ortholog),
                          Single_copy_orthologs = as.numeric(Single_copy_orthologs)
)
str(genome_stats)
#write.table(genome_stats, "../01_funannotate/funannotate_compare_sporothrix_20220915/stats/genome.stats.summary_wide.txt", sep='\t', quote = F, row.names = T)
genome_stats_matrix <- as.matrix(genome_stats)

assembly_size_matrix <- genome_stats_matrix[,1]
no_proteins_matrix <- genome_stats_matrix[,8]

secmets<-read.delim("antiSMASH/antismash.txt") %>% dplyr::select(genome, secondary_metabolite_type) %>%
  dplyr::group_by(genome) %>%
  dplyr::count(secondary_metabolite_type) %>%
  tidyr::spread(secondary_metabolite_type, n, fill=0) %>%
  base::as.data.frame()
secmets <- left_join(genome_df,secmets) %>%
  select(genome, NRPS,T1PKS,T3PKS,terpene,"fungal-RiPP-like",other)
row.names(secmets) = secmets$genome
secmets <- secmets %>% dplyr::select(-genome)
secmets_matrix <- as.matrix(secmets)

CAZyme <- read.csv("funannotate/funannotate_compare_20230729/cazy/CAZyme.summary.results.csv")
row.names(CAZyme) = CAZyme$X
CAZyme <- CAZyme %>% dplyr::select(-X)
CAZyme <- as.data.frame(t(CAZyme)) #transpose
CAZyme <- CAZyme %>% add_column(genome = rownames(CAZyme), .before = "AA") #make new column with the row names
CAZyme <- left_join(genome_df,CAZyme)
row.names(CAZyme) = CAZyme$genome

CAZyme <- CAZyme %>% dplyr::select(-genome)
CAZyme_matrix <- as.matrix(CAZyme)

merops<- read.csv("funannotate/funannotate_compare_20230729/merops/MEROPS.summary.results.csv")
row.names(merops) = merops$X
merops <- merops %>% dplyr::select(-X)
merops <- as.data.frame(t(merops))
merops <- merops %>% add_column(genome = rownames(merops), .before = "A") 
merops <- left_join(genome_df,merops)
row.names(merops) = merops$genome
merops <- merops %>% dplyr::select(-genome)
merops_matrix <- as.matrix(merops)


# make heatmaps and annotations

annotations <- rowAnnotation(
  "compl" = anno_points((completeness_matrix),
                        width = unit(2, "cm"),
                        ylim = c(0, 100),
                        pch = 1:2,
                        axis_param = list(
                          at = c(0, 50, 100), 
                          labels = c("0", "50", "100")),
                        #gp = gpar(col = 2:3),
                        border = FALSE),
  "bp" = anno_barplot((assembly_size_matrix),
                      width = unit(1, "cm"),
                      baseline = "min",
                      ylim = c(18000000, 44000000),
                      axis_param = list(
                        at = c(20000000, 42000000), 
                        labels = c("20", "42")),
                      border = FALSE,
                      gp = gpar(fill = "#2E3033", col = "#2E3033", lineend = "round")),
  "proteins" = anno_barplot((no_proteins_matrix),
                            width = unit(1, "cm"),
                            baseline = "min",
                            ylim = c(6000, 11000),
                            axis_param = list(
                              at = c(7000, 10000), 
                              labels = c("7k", "10k")),
                            border = FALSE,
                            gp = gpar(fill = "#2E3033", col = "#2E3033", lineend = "round")),
  gap = unit(3, "mm"),
  show_annotation_name = FALSE)

library(circlize)
col_fun = colorRamp2(c(0, max(unlist(CAZyme_matrix))), c("white", "darkorchid"))
col_fun(seq(0, 3))
cazyme_hm <- Heatmap(CAZyme_matrix,
                     name = "CAZymes",
                     column_title = "CAZyme\nclasses",
                     column_title_gp = gpar(fontsize = 12),
                     col = col_fun,
                     cluster_rows = ML_dend,
                     row_dend_width = unit(2, "cm"),
                     column_names_gp = grid::gpar(fontsize = 8),
                     row_names_side = "left",
                     row_names_gp = grid::gpar(fontsize = 8),
                     width = unit(3, "cm"),
                     height = unit(8, "cm"),
                     left_annotation = annotations,
                     #left_annotation = rowAnnotation(
                     #pathogen = genome_attributes_matrix[,"mammal_pathogen"],
                     #soil = genome_attributes_matrix[,"soil"],
                     #arthropod = genome_attributes_matrix[,"arthropod"],
                     #wood = genome_attributes_matrix[,"wood"],
                     #infructescence = genome_attributes_matrix[,"infructescence"],
                     #macrofungus = genome_attributes_matrix[,"macrofungus"],
                     #col = list(arthropod = c("1" = "black",
                     #"NA" = "azure"),
                     # pathogen = c("1" = "black",
                     #"NA" = "white"),
                     #soil = c("1" = "black",
                     #"NA" = "white"),
                     #wood = c("1" = "black",
                     #"NA" = "white"),                                  
                     #infructescence = c("1" = "black",
                     #"NA" = "white"),                                  
                     #macrofungus = c("1" = "black",
                     #"NA" = "white")),
                     #show_legend = FALSE),
                     show_heatmap_legend = FALSE,
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       grid.text(CAZyme_matrix[i, j], x, y, gp = gpar(fontsize = 6))})
cazyme_hm

col_fun = colorRamp2(c(0, max(unlist(secmets_matrix))), c("white", "deepskyblue4"))
secmets_hm <- Heatmap(secmets_matrix,
                      name = "BGCs",
                      col = col_fun,
                      column_title = "Biosynthetic gene\nclusters",
                      column_title_gp = gpar(fontsize = 12),
                      cluster_rows = ML_dend,
                      row_dend_width = unit(2, "cm"),
                      column_names_gp = grid::gpar(fontsize = 8),
                      row_names_side = "left",
                      row_names_gp = grid::gpar(fontsize = 8),
                      width = unit(3.5, "cm"),
                      height = unit(8, "cm"),
                      #right_annotation = annotations,
                      show_heatmap_legend = FALSE,
                      cell_fun = function(j, i, x, y, width, height, fill) {
                        grid.text(secmets_matrix[i, j], x, y, gp = gpar(fontsize = 6))})
secmets_hm

col_fun = colorRamp2(c(0, max(unlist(merops_matrix))), c("white", "aquamarine4"))
merops_hm <- Heatmap(merops_matrix,
                     name = "proteases",
                     column_title = "Protease\nclasses",
                     column_title_gp = gpar(fontsize = 12),
                     col = col_fun,
                     cluster_rows = ML_dend,
                     row_dend_width = unit(2, "cm"),
                     column_names_gp = grid::gpar(fontsize = 8),
                     row_names_side = "left",
                     row_names_gp = grid::gpar(fontsize = 8),
                     width = unit(4.5, "cm"),
                     height = unit(8, "cm"),
                     #right_annotation = annotations,
                     show_heatmap_legend = FALSE,
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       grid.text(merops_matrix[i, j], x, y, gp = gpar(fontsize = 6))})
merops_hm

library(circlize)
col_fun = colorRamp2(c(0, 1), c("white", "grey"))
col_fun(seq(0, 3))

attributes_hm <- Heatmap(genome_attributes_lifestyle_matrix,
                         name = "Attributes",
                         column_title = "Association\n",
                         column_title_gp = gpar(fontsize = 12),
                         col = col_fun,
                         cluster_rows = ML_dend,
                         row_dend_width = unit(2, "cm"),
                         column_names_gp = grid::gpar(fontsize = 8),
                         row_names_side = "left",
                         row_names_gp = grid::gpar(fontsize = 8),
                         width = unit(3, "cm"),
                         height = unit(8, "cm"),
                         show_heatmap_legend = FALSE)
attributes_hm 


attributes_hm + cazyme_hm + secmets_hm + merops_hm
decorate_annotation("compl", { 
  grid.text("Assembly", y = unit(1, "npc") + unit(18, "mm"), just = "bottom", rot = 0, gp = gpar(fontsize = 12)) 
})
decorate_annotation("compl", { 
  grid.text("completeness", y = unit(1, "npc") + unit(13, "mm"), just = "bottom", rot = 0, gp = gpar(fontsize = 12)) 
})

decorate_annotation("bp", { 
  grid.text("Mbp", y = unit(1, "npc") + unit(18, "mm"), just = "bottom", rot = 0, gp = gpar(fontsize = 12)) 
})

decorate_annotation("proteins", { 
  grid.text("Predicted", y = unit(1, "npc") + unit(18, "mm"), just = "bottom", rot = 0, gp = gpar(fontsize = 12)) 
})
decorate_annotation("proteins", { 
  grid.text("proteins", y = unit(1, "npc") + unit(13, "mm"), just = "bottom", rot = 0, gp = gpar(fontsize = 12)) 
})

#write to pdf
pdf(file="~/Documents/2022_04_eukaryote_annotation/03_heatmap_annotations/absolute_values_sporothrix.pdf", width=12, height=8)
attributes_hm + cazyme_hm + secmets_hm + merops_hm
decorate_annotation("compl", { 
  grid.text("Assembly", y = unit(1, "npc") + unit(18, "mm"), just = "bottom", rot = 0, gp = gpar(fontsize = 12)) 
})
decorate_annotation("compl", { 
  grid.text("completeness", y = unit(1, "npc") + unit(13, "mm"), just = "bottom", rot = 0, gp = gpar(fontsize = 12)) 
})

decorate_annotation("bp", { 
  grid.text("Mbp", y = unit(1, "npc") + unit(18, "mm"), just = "bottom", rot = 0, gp = gpar(fontsize = 12)) 
})

decorate_annotation("proteins", { 
  grid.text("Predicted", y = unit(1, "npc") + unit(18, "mm"), just = "bottom", rot = 0, gp = gpar(fontsize = 12)) 
})
decorate_annotation("proteins", { 
  grid.text("proteins", y = unit(1, "npc") + unit(13, "mm"), just = "bottom", rot = 0, gp = gpar(fontsize = 12)) 
})
dev.off()