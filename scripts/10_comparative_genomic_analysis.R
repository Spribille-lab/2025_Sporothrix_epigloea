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

# make a heatmap using complexheatmap that compares Sporothrix species
setwd("~/Documents/2023_02_Sporothrix/10_annotations/")

ML_tree <- ape::read.tree("~/Documents/2023_02_Sporothrix/09_orthology/OrthoFinder/Results_Aug03/Species_Tree/iqtree2/sporothrix_placement.rerooted.contree")
ML_tree$tip.label
ML_tree <- drop.tip(ML_tree, 25)
str(ML_tree)

ML_dend <- ReadDendrogram(textConnection(write.tree(ML_tree)))
str(ML_dend)

#ML_dend <- DECIPHER::ReadDendrogram("~/Documents/2022_04_eukaryote_annotation/05_genomic_analysis_Sporothrix/ML_tree2.txt",keepRoot = T)
order.dendrogram(ML_dend)
order.dendrogram(ML_dend) <- 1:nleaves(ML_dend)
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

unannotated_genes <- read.delim("unannotated_genes.txt")

genome_stats <-read.csv("funannotate/funannotate_compare_20230809/stats/genome.stats.summary.csv")
row.names(genome_stats) = genome_stats$X
genome_stats <- genome_stats %>% dplyr::select(-X)
genome_stats <- as.data.frame(t(genome_stats)) #transpose
rownames(genome_stats)
genome_stats <- genome_stats %>% add_column(genome = rownames(genome_stats), .before = "isolate") #make new column with the row names
genome_stats <- genome_df %>% left_join(genome_stats) %>% left_join(unannotated_genes) %>% select(-genes)
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
                          Single_copy_orthologs = as.numeric(Single_copy_orthologs),
                          Percent_masked_repeats = as.numeric(Percent_masked_repeats),
                          unannotated = as.numeric(unannotated)
)
str(genome_stats)
#write.table(genome_stats, "../01_funannotate/funannotate_compare_sporothrix_20220915/stats/genome.stats.summary_wide.txt", sep='\t', quote = F, row.names = T)
genome_stats_matrix <- as.matrix(genome_stats)




assembly_size_matrix <- genome_stats_matrix[,1]
no_proteins_matrix <- genome_stats_matrix[,8]
repeats_matrix <- genome_stats_matrix[,13]
N50_matrix <- genome_stats_matrix[,5]
scaffolds_matrix <- genome_stats_matrix[,4]
unannotated_matrix <- genome_stats_matrix[,c(14,7)]
proportion_annotated_matrix <- t(apply(unannotated_matrix, 1, function(x) x/sum(x)))

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

CAZyme <- read.csv("funannotate/funannotate_compare_20230929/cazy/CAZyme.summary.results.csv")
row.names(CAZyme) = CAZyme$X
CAZyme <- CAZyme %>% dplyr::select(-X)
CAZyme <- as.data.frame(t(CAZyme)) #transpose
CAZyme <- CAZyme %>% add_column(genome = rownames(CAZyme), .before = "AA") #make new column with the row names
CAZyme <- left_join(genome_df,CAZyme)
row.names(CAZyme) = CAZyme$genome

CAZyme <- CAZyme %>% dplyr::select(-genome)
CAZyme_matrix <- as.matrix(CAZyme)

merops<- read.csv("funannotate/funannotate_compare_20230809/merops/MEROPS.summary.results.csv")
row.names(merops) = merops$X
merops <- merops %>% dplyr::select(-X)
merops <- as.data.frame(t(merops))
merops <- merops %>% add_column(genome = rownames(merops), .before = "A") 
merops <- left_join(genome_df,merops)
row.names(merops) = merops$genome
merops <- merops %>% dplyr::select(-genome)
merops_matrix <- as.matrix(merops)

#YR<-read.delim("~/Documents/2023_02_Sporothrix/12_starfish/tyrCount.txt")
#YR <- left_join(genome_df,YR)
#row.names(YR) = YR$genome
#YR <- YR %>% dplyr::select(-genome)
#YR_matrix <- as.matrix(YR)

# make heatmaps and annotations


row_labels = genome
row_labels
names(row_labels) = c("S. euskadiensis VPRI43754",
                      "S. pseudoabietina VPRI43531",
                      "S. variecibatus CBS 121961",
                      "S. protearum CBS 116654",
                      "S. humicola CBS 118129",
                      "S. pallida CBS 131.56",
                      "S. mexicana CBS 120341",
                      "S. brasiliensis 5110",
                      "S. schenckii 1099",
                      "S. globosa CBS 120340",
                      "S. luriei CBS 937.72",
                      "S. phasma CBS 119721",
                      "S. dimorphospora CBS 553.74",
                      "S. inflata CBS 239.68",
                      "S. bragantina CBS 474.91",
                      "S. curviconia CBS 959.73",
                      "S. thermara CBS 139747MAG",
                      "S. epigloea CBS 573.63",
                      "S. epigloea TF4163MAG",
                      "S. epigloea CBS 119000",
                      "S. eucalyptigena CBS 139899",
                      "S. eucalyptigena CBS 140593",
                      "S. nigrograna VPRI43755",
                      "O. novo-ulmi H327",
                      "O. ips VPRI43529",
                      "S. brunneoviolacea CBS 124561",
                      "O. fasciatum VPRI43845",
                      "S. insectorum RCEF264",
                      "L. lundbergii CBS 138716")
names(row_labels)




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
  "N50" = anno_text((N50_matrix),
                         location = 0.5,
                         just = "center",
                         #width = unit(2, "cm"),
                         #border = FALSE,
                         gp = gpar(fontsize = 7, col = "black", border = "black", fill = "aliceblue"),
                         width = max_text_width(N50_matrix)*1.05),
  "contigs" = anno_text((scaffolds_matrix),
                        location = 0.5,
                        just = "center",
                        #width = unit(2, "cm"),
                        #border = FALSE,
                        gp = gpar(fontsize = 7, col = "black", border = "black", fill = "cornsilk1"),
                        width = max_text_width(scaffolds_matrix)*1.05),
  "bp" = anno_barplot((assembly_size_matrix),
                      width = unit(1, "cm"),
                      baseline = c(18000000),
                      ylim = c(18000000, 44000000),
                      axis_param = list(
                        at = c(20000000, 42000000), 
                        labels = c("20", "42")),
                      border = FALSE,
                      gp = gpar(fill = "#2E3033", col = "#2E3033", lineend = "round")),
  "proteins" = anno_barplot((no_proteins_matrix),
                            width = unit(1, "cm"),
                            baseline = c(6000),
                            ylim = c(6000, 11000),
                            axis_param = list(
                              at = c(7000, 10000), 
                              labels = c("7k", "10k")),
                            border = FALSE,
                            gp = gpar(fill = "#2E3033", col = "#2E3033", lineend = "round")),
  "repeats" = anno_barplot((repeats_matrix),
                      width = unit(1, "cm"),
                      baseline = c(0),
                      ylim = c(0, 18),
                      axis_param = list(
                        at = c(0, 15), 
                        labels = c("0", "15")),
                      border = FALSE,
                      gp = gpar(fill = "#2E3033", col = "#2E3033", lineend = "round")),
  "contigs" = anno_barplot((scaffolds_matrix),
                           width = unit(2, "cm"),
                           baseline = c(0),
                           #ylim = c(0, 18),
                          # axis_param = list(
                          #   at = c(0, 15), 
                           #  labels = c("0", "15")),
                           border = FALSE,
                           gp = gpar(fill = "#2E3033", col = "#2E3033", lineend = "round")),
  gap = unit(4, "mm"),
  show_annotation_name = FALSE)

library(circlize)
col_fun = colorRamp2(c(0, 1), c("white", "gray44"))
col_fun(seq(0, 3))

attributes_hm <- Heatmap(genome_attributes_lifestyle_matrix,
                         name = "Attributes",
                         column_title = "Association\n",
                         column_title_gp = gpar(fontsize = 10),
                         col = col_fun,
                         cluster_rows = ML_dend,
                         row_dend_width = unit(2, "cm"),
                         column_names_gp = grid::gpar(fontsize = 8),
                         row_names_side = "left",
                         row_labels = names(row_labels),
                         row_names_gp = grid::gpar(fontsize = 8),
                         width = unit(3, "cm"),
                         height = unit(9, "cm"),
                         show_heatmap_legend = FALSE,
                         show_column_dend = FALSE,
                         right_annotation = annotations)

attributes_hm 

#write to pdf
pdf(file="~/Documents/2023_02_Sporothrix/results/figures/sporothrix_assembly_details.pdf", width=10, height=8)
attributes_hm
decorate_annotation("compl", { 
  grid.text("Assembly", y = unit(1, "npc") + unit(12, "mm"), just = "bottom", rot = 0, gp = gpar(fontsize = 10)) 
})
decorate_annotation("compl", { 
  grid.text("completeness", y = unit(1, "npc") + unit(7, "mm"), just = "bottom", rot = 0, gp = gpar(fontsize = 10)) 
})

decorate_annotation("bp", { 
  grid.text("Mbp", y = unit(1, "npc") + unit(7, "mm"), just = "bottom", rot = 0, gp = gpar(fontsize = 10)) 
})

decorate_annotation("proteins", { 
  grid.text("Predicted", y = unit(1, "npc") + unit(12, "mm"), just = "bottom", rot = 0, gp = gpar(fontsize = 10)) 
})
decorate_annotation("proteins", { 
  grid.text("proteins", y = unit(1, "npc") + unit(7, "mm"), just = "bottom", rot = 0, gp = gpar(fontsize = 10)) 
})
decorate_annotation("repeats", { 
  grid.text("repeats", y = unit(1, "npc") + unit(7, "mm"), just = "bottom", rot = 0, gp = gpar(fontsize = 10)) 
})
decorate_annotation("N50", { 
  grid.text("N50", y = unit(1, "npc") + unit(7, "mm"), just = "bottom", rot = 0, gp = gpar(fontsize = 10)) 
})
decorate_annotation("contigs", { 
  grid.text("Contigs", y = unit(1, "npc") + unit(7, "mm"), just = "bottom", rot = 0, gp = gpar(fontsize = 10)) 
})




dev.off()



class_annotations <- rowAnnotation(
  "unannotated" = anno_barplot((unannotated_matrix[,1]),
                            width = unit(1.5, "cm"),
                            #baseline = c(6000),
                            #ylim = c(6000, 11000),
                            axis_param = list(
                              at = c(1000, 2000), 
                              labels = c("1k", "2k")),
                            border = FALSE,
                            gp = gpar(fill = "#2E3033", col = "#2E3033", lineend = "round")),
  "proportion" = anno_barplot(proportion_annotated_matrix,
                              width = unit(2, "cm"),
                              gp = gpar(fill = c("#2E3033","lightgrey"), col = "#2E3033"),
                              #bar_width = 2,
                              #gap = unit(4, "mm")
                              border = FALSE),
  gap = unit(4, "mm"),
  show_annotation_name = FALSE)


library(circlize)
col_fun = colorRamp2(c(0, max(unlist(CAZyme_matrix))), c("white", "darkorchid"))
col_fun(seq(0, 3))
cazyme_hm <- Heatmap(CAZyme_matrix,
                     name = "CAZymes",
                     column_title = "CAZyme\nclasses",
                     column_title_gp = gpar(fontsize = 10),
                     col = col_fun,
                     cluster_rows = ML_dend,
                     row_dend_width = unit(2, "cm"),
                     column_names_gp = grid::gpar(fontsize = 8),
                     row_names_side = "left",
                     row_labels = names(row_labels),
                     row_names_gp = grid::gpar(fontsize = 8),
                     width = unit(3, "cm"),
                     height = unit(8, "cm"),
                     show_column_dend = FALSE,
                     #right_annotation = class_annotations,
                     show_heatmap_legend = FALSE,
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       grid.text(CAZyme_matrix[i, j], x, y, gp = gpar(fontsize = 6))})
cazyme_hm

col_fun = colorRamp2(c(0, max(unlist(secmets_matrix))), c("white", "aquamarine4"))
secmets_hm <- Heatmap(secmets_matrix,
                      name = "BGCs",
                      col = col_fun,
                      column_title = "Biosynthetic gene\nclusters",
                      column_title_gp = gpar(fontsize = 10),
                      cluster_rows = ML_dend,
                      row_dend_width = unit(2, "cm"),
                      column_names_gp = grid::gpar(fontsize = 8),
                      row_names_side = "left",
                      row_labels = names(row_labels),
                      row_names_gp = grid::gpar(fontsize = 8),
                      width = unit(3.5, "cm"),
                      height = unit(8, "cm"),
                      right_annotation = class_annotations,
                      show_heatmap_legend = FALSE,
                      show_column_dend = FALSE,
                      cell_fun = function(j, i, x, y, width, height, fill) {
                        grid.text(secmets_matrix[i, j], x, y, gp = gpar(fontsize = 6))})
secmets_hm

col_fun = colorRamp2(c(0, max(unlist(merops_matrix))), c("white", "deepskyblue4"))
merops_hm <- Heatmap(merops_matrix,
                     name = "proteases",
                     column_title = "Protease\nclasses",
                     column_title_gp = gpar(fontsize = 10),
                     col = col_fun,
                     cluster_rows = ML_dend,
                     row_dend_width = unit(2, "cm"),
                     column_names_gp = grid::gpar(fontsize = 8),
                     row_names_side = "left",
                     row_labels = names(row_labels),
                     row_names_gp = grid::gpar(fontsize = 8),
                     width = unit(4.5, "cm"),
                     height = unit(8, "cm"),
                     show_heatmap_legend = FALSE,
                     show_column_dend = FALSE,
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       grid.text(merops_matrix[i, j], x, y, gp = gpar(fontsize = 6))})
merops_hm

#write to pdf
pdf(file="~/Documents/2023_02_Sporothrix/results/figures/sporothrix_class_inventory.pdf", width=12, height=6)
cazyme_hm + merops_hm + secmets_hm
decorate_annotation("unannotated", { 
  grid.text("counts", y = unit(1, "npc") + unit(2, "mm"), just = "bottom", rot = 0, gp = gpar(fontsize = 10)) 
})
decorate_annotation("unannotated", { 
  grid.text("Unannotated genes", y = unit(1, "npc") + unit(10, "mm"), just = "bottom", rot = 0, gp = gpar(fontsize = 10)) 
})
decorate_annotation("proportion", { 
  grid.text("proportion\nof total", y = unit(1, "npc") + unit(3, "mm"), just = "bottom", rot = 0, gp = gpar(fontsize = 10)) 
})


dev.off()



















##########calculate the increase and decrease for epigloea##########

group = c(rep("other_sporothrix",17), rep("epigloea",3), rep("other_sporothrix",3), rep("outgroup", 6))

CAZyme <- CAZyme %>% mutate(group = group) %>% filter(!group %in% c("outgroup"))
CAZyme_tidy <- CAZyme %>% gather("AA", "CBM", "CE", "GH", "GT", "PL", key = "CAZyme_class", value = "n")

wilcox.test(AA ~ group, data = CAZyme)
wilcox.test(CBM ~ group, data = CAZyme)
wilcox.test(CE ~ group, data = CAZyme)
wilcox.test(GH ~ group, data = CAZyme)
wilcox.test(GT ~ group, data = CAZyme)
wilcox.test(PL ~ group, data = CAZyme)

CAZyme_stats <- CAZyme_tidy %>%
  group_by(CAZyme_class) %>%
  group_by(group, .add = TRUE) %>%
  summarize(mean = mean(n)) %>%
  spread(group, mean) %>%
  ungroup()

write.table(CAZyme_stats, "funannotate/increase_decrease_stats/CAZyme.sporothrix.percent.txt", sep='\t', quote = F, row.names = F)


secmets <- secmets %>% mutate(group = group) %>% filter(!group %in% c("outgroup"))
colnames(secmets) = c("NRPS", "T1PKS", "T3PKS", "terpene", "RiPP", "other", "group")
secmets_tidy <- secmets %>% gather("NRPS", "T1PKS", "T3PKS", "terpene", "RiPP", "other", key = "BGC_class", value = "n")
secmets_stats <- secmets_tidy %>%
  group_by(BGC_class) %>%
  group_by(group, .add = TRUE) %>%
  summarize(mean = mean(n)) %>%
  spread(group, mean) %>%
  ungroup()

write.table(secmets_stats, "funannotate/increase_decrease_stats/secmets.sporothrix.percent.txt", sep='\t', quote = F, row.names = F)

wilcox.test(NRPS ~ group, data = secmets)
wilcox.test(T1PKS ~ group, data = secmets)
wilcox.test(T3PKS ~ group, data = secmets)
wilcox.test(terpene ~ group, data = secmets)
wilcox.test(RiPP ~ group, data = secmets)
wilcox.test(other ~ group, data = secmets)

merops <- merops %>% mutate(group = group) %>% filter(!group %in% c("outgroup"))
merops_tidy <- merops %>% gather("A", "C", "G", "M", "N", "P", "S", "T", "I", key = "protease_class", value = "n")
merops_stats <- merops_tidy %>%
  group_by(protease_class) %>%
  group_by(group, .add = TRUE) %>%
  summarize(mean = mean(n)) %>%
  spread(group, mean)

write.table(merops_stats, "funannotate/increase_decrease_stats/merops.sporothrix.percent.txt", sep='\t', quote = F, row.names = F)

wilcox.test(A ~ group, data = merops)
wilcox.test(C ~ group, data = merops)
wilcox.test(G ~ group, data = merops)
wilcox.test(M ~ group, data = merops)
wilcox.test(N ~ group, data = merops)
wilcox.test(P ~ group, data = merops)
wilcox.test(S ~ group, data = merops)
wilcox.test(T ~ group, data = merops)
wilcox.test(I ~ group, data = merops)

