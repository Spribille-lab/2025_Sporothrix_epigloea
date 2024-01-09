library(tidyverse)
library(data.table)
library(dplyr)
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


############## CAZyme data ###############
setwd("~/Documents/2023_02_Sporothrix/10_annotations/")

# secreted CAZymes

file_list = list.files(path="annotation_tables_merged/")

cazymes_secreted <- data.frame()
for (i in 1:length(file_list)){
  temp_data <- read.delim(paste("annotation_tables_merged/",file_list[i], sep = ""), stringsAsFactors = F)
  temp_data <- temp_data %>%
    select(Gene_ID,Localizations_DeepLoc2,CAZyme) %>%
    add_column(genome = file_list[i], .before = "Gene_ID") %>%
    dplyr::mutate(genome = str_replace(genome, "(\\S*)(.annotations.txt)", "\\1")) %>%
    dplyr::mutate(organism = str_replace(genome, "(\\S*)(_)(\\S*)(_)(\\S*_)(\\S*)(_)(\\S*_)(\\S*)", "\\1\\2\\3\\4\\6")) %>%
    dplyr::relocate(organism) %>%
    separate_rows(CAZyme, sep = ";") %>%
    filter(Localizations_DeepLoc2 == "Extracellular")
  temp_data <- temp_data[!is.na(temp_data$CAZyme),]
  cazymes_secreted <- rbindlist(list(cazymes_secreted, temp_data), use.names = T)
}

head(cazymes_secreted)

#######################################################################

# count the number of hits per genome
cazyme_count <- cazymes_secreted %>%
  dplyr::group_by(organism) %>%
  dplyr::count(CAZyme)
colnames(cazyme_count) = c("organism", "CAZyme_family", "n")

# make a wide dataframe of the counts
cazyme_count_wide <- cazyme_count %>%
  tidyr::spread(organism, n, fill=0) %>%
  base::as.data.frame() %>%
  dplyr::arrange(CAZyme_family)

#make a list of the CAZyme families and their putative functions

CAZyme_activity <- read.delim("~/Documents/2022_04_eukaryote_annotation/06_dbcan/CAZyme_activity.txt", header = T, na.strings='')
CAZyme_activity <- CAZyme_activity %>% select(General_Activity, CAZyme_families, origin)
colnames(CAZyme_activity) = c("General_Activity", "CAZyme_family", "origin")


cazyme_count_wide <- cazyme_count_wide %>%
  dplyr::left_join(CAZyme_activity, by = "CAZyme_family") %>%
  dplyr::arrange(General_Activity)

write.table(cazyme_count_wide, "tables/sporothrix_secreted_CAZymes_counts.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)  


# make matrix for full heatmap
cazyme_matrix <- cazyme_count_wide %>%
  select(-General_Activity, -origin)
row.names(cazyme_matrix) = cazyme_matrix$CAZyme_family
cazyme_matrix <- cazyme_matrix %>% dplyr::select(-CAZyme_family)
cazyme_matrix <- as.matrix(cazyme_matrix)

cazyme_matrix_t <- cazyme_count_wide %>%
  select(-General_Activity, -origin) %>%
  arrange(CAZyme_family)
row.names(cazyme_matrix_t) = cazyme_matrix_t$CAZyme_family
cazyme_matrix_t <- cazyme_matrix_t  %>% select(-CAZyme_family)
cazyme_matrix_t <- as.data.frame(t(cazyme_matrix_t))
cazyme_matrix_t <- cazyme_matrix_t %>% rownames_to_column(var = "genome")
cazyme_matrix_t <- left_join(genome_df, cazyme_matrix_t)
row.names(cazyme_matrix_t) = cazyme_matrix_t$genome
cazyme_matrix_t <- cazyme_matrix_t %>% dplyr::select(-genome)
cazyme_matrix_t <- as.matrix(cazyme_matrix_t)


CAZyme_activity_matrix <- cazyme_count_wide %>% select(CAZyme_family, General_Activity, origin)
row.names(CAZyme_activity_matrix) = CAZyme_activity_matrix$CAZyme_family
CAZyme_activity_matrix <- CAZyme_activity_matrix %>%
  dplyr::select(-CAZyme_family) %>%
  as.matrix()

################ vertical heatmap #############

library(circlize)
col_fun = colorRamp2(c(0, max(unlist(cazyme_matrix))), c("white", "darkorchid"))
max(unlist(cazyme_matrix))

family_annotation <- rowAnnotation(
  "General_activity" = anno_text(CAZyme_activity_matrix[,1], gp = gpar(fontsize = 6)),
  "origin" = anno_text(CAZyme_activity_matrix[,2], gp = gpar(fontsize = 6)))

cazyme_heatmap <- Heatmap(cazyme_matrix,
                          name = "n",
                          row_order = order(CAZyme_activity_matrix[, 2]),
                          column_names_gp = grid::gpar(fontsize = 7),
                          row_names_gp = grid::gpar(fontsize = 6),
                          row_names_side = "left",
                          #cluster_rows = FALSE,
                          col = col_fun,
                          width = unit(9, "cm"),
                          height = unit(30, "cm"),
                          right_annotation = family_annotation,
                          show_heatmap_legend = FALSE,
                          cell_fun = function(j, i, x, y, width, height, fill) {
                            grid.text(cazyme_matrix[i, j], x, y, gp = gpar(fontsize = 6))}
)
draw(cazyme_heatmap)

pdf(file="figures/sporothrix_cazyme_heatmap_complexheatmap.pdf", width=10, height=15)

draw(cazyme_heatmap)
dev.off()

# make a list of all the secreted CAZymes that I need to annotate

secreted_CAZymes_sporothrix <- cazymes_secreted %>%
  distinct(CAZyme) %>%
  pull(CAZyme)
write.table(secreted_CAZymes_sporothrix, "tables/sporothrix_secreted_CAZymes.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)  

################ horizontal heatmap #################
library(circlize)
col_fun = colorRamp2(c(0, max(unlist(cazyme_matrix_t))), c("white", "darkorchid"))
max(unlist(cazyme_matrix_t))

wide_heatmap <- Heatmap(cazyme_matrix_t,
                        col = col_fun,
                        name = "n",
                        column_names_gp = grid::gpar(fontsize = 7),
                        row_names_gp = grid::gpar(fontsize = 7),
                        row_names_side = "left",
                        cluster_rows = ML_dend,
                        cluster_columns = FALSE,
                        width = unit(30, "cm"),
                        height = unit(10, "cm"),
                        show_heatmap_legend = FALSE,
                          cell_fun = function(j, i, x, y, width, height, fill) {
                            grid.text(cazyme_matrix_t[i, j], x, y, gp = gpar(fontsize = 6))}
)

#row_order = order(CAZyme_activity_matrix[, 2]),
draw(wide_heatmap)

pdf(file="figures/sporothrix_cazyme_heatmap_wide_complexheatmap.pdf", width=18, height=8)

draw(wide_heatmap)
dev.off()
############# protease data ############

setwd("~/Documents/2023_02_Sporothrix/10_annotations/")

# secreted CAZymes

file_list = list.files(path="annotation_tables_merged/")

proteases_secreted <- data.frame()
for (i in 1:length(file_list)){
  temp_data <- read.delim(paste("annotation_tables_merged/",file_list[i], sep = ""), stringsAsFactors = F)
  temp_data <- temp_data %>%
    select(Gene_ID,Localizations_DeepLoc2,Protease) %>%
    add_column(genome = file_list[i], .before = "Gene_ID") %>%
    dplyr::mutate(genome = str_replace(genome, "(\\S*)(.annotations.txt)", "\\1")) %>%
    dplyr::mutate(organism = str_replace(genome, "(\\S*)(_)(\\S*)(_)(\\S*_)(\\S*)(_)(\\S*_)(\\S*)", "\\1\\2\\3\\4\\6")) %>%
    dplyr::relocate(organism) %>%
    separate_rows(Protease, sep = ";") %>%
    filter(Localizations_DeepLoc2 == "Extracellular")
  temp_data <- temp_data[!is.na(temp_data$Protease),]
  proteases_secreted <- rbindlist(list(proteases_secreted, temp_data), use.names = T)
}

head(proteases_secreted)

# count the number of hits per genome
protease_count <- proteases_secreted %>%
  dplyr::group_by(organism) %>%
  dplyr::count(Protease)
colnames(protease_count) = c("organism", "protease_family", "n")

# make a wide dataframe of the counts
protease_count_wide <- protease_count %>%
  tidyr::spread(organism, n, fill=0) %>%
  base::as.data.frame() %>%
  dplyr::arrange(protease_family)


write.table(protease_count_wide, "tables/sporothrix_secreted_proteases_counts.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)  


# make matrix for full heatmap
protease_matrix <- protease_count_wide
row.names(protease_matrix) = protease_matrix$protease_family
protease_matrix <- protease_matrix %>% dplyr::select(-protease_family)
protease_matrix <- as.matrix(protease_matrix)

protease_matrix_t <- protease_count_wide %>%
  arrange(protease_family)
row.names(protease_matrix_t) = protease_matrix_t$protease_family
protease_matrix_t <- protease_matrix_t %>% select(-protease_family)
protease_matrix_t <- as.data.frame(t(protease_matrix_t))
protease_matrix_t <- protease_matrix_t %>% rownames_to_column(var = "genome")
protease_matrix_t <- left_join(genome_df, protease_matrix_t)
row.names(protease_matrix_t) = protease_matrix_t$genome
protease_matrix_t <- protease_matrix_t %>% dplyr::select(-genome)
protease_matrix_t <- as.matrix(protease_matrix_t)

################ vertical heatmap #############

library(circlize)
col_fun = colorRamp2(c(0, max(unlist(protease_matrix))), c("white", "aquamarine4"))
max(unlist(protease_matrix))

protease_heatmap <- Heatmap(protease_matrix,
                          name = "n",
                          #row_order = order(protease_activity_matrix[, 2]),
                          column_names_gp = grid::gpar(fontsize = 7),
                          row_names_gp = grid::gpar(fontsize = 6),
                          row_names_side = "left",
                          #cluster_rows = FALSE,
                          col = col_fun,
                          width = unit(9, "cm"),
                          height = unit(10, "cm"),
                          show_heatmap_legend = FALSE,
                          cell_fun = function(j, i, x, y, width, height, fill) {
                            grid.text(protease_matrix[i, j], x, y, gp = gpar(fontsize = 6))}
)
draw(protease_heatmap)

pdf(file="figures/sporothrix_protease_heatmap_complexheatmap.pdf", width=10, height=15)

draw(protease_heatmap)
dev.off()

################ horizontal heatmap #################
library(circlize)
col_fun = colorRamp2(c(0, max(unlist(protease_matrix_t))), c("white", "aquamarine4"))
max(unlist(protease_matrix))

wide_heatmap <- Heatmap(protease_matrix_t,
                        col = col_fun,
                        name = "n",
                        column_names_gp = grid::gpar(fontsize = 7),
                        row_names_gp = grid::gpar(fontsize = 7),
                        row_names_side = "left",
                        cluster_rows = ML_dend,
                        #cluster_columns = FALSE,
                        width = unit(10, "cm"),
                        height = unit(10, "cm"),
                        show_heatmap_legend = FALSE,
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(protease_matrix_t[i, j], x, y, gp = gpar(fontsize = 6))}
)

#row_order = order(protease_activity_matrix[, 2]),
draw(wide_heatmap)

pdf(file="figures/sporothrix_protease_heatmap_wide_complexheatmap.pdf", width=12, height=8)

draw(wide_heatmap)
dev.off()
