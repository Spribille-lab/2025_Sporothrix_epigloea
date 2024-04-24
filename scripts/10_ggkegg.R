#BiocManager::install("ggkegg")
library(ggkegg)
library(dplyr)
library(tidygraph)
library(ggfx)
library(ggraph)
library(igraph)
library(clusterProfiler)
library(kableExtra)

# Assessing module completeness across multiple microbial genomes
setwd("~/Documents/2023_02_Sporothrix/10_annotations/kegg/")

############## get tree #############
library("ape")
library('dendextend')
library("DECIPHER")

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
######################


## get module files from the kegg API

# cd /Users/carmenallen/Documents/2023_02_Sporothrix/10_annotations/kegg
# wget https://rest.kegg.jp/list/module
# cd mf
# while IFS=$'\t' read -r v1 v2; do wget https://rest.kegg.jp/get/"$v1"; sleep 1; done < ../module
# 

## Load pre-computed module files obtained from the kegg API
mf <- list.files("mf")
annos <- list()
genomes <- read.delim("../genomes.txt")
genomes <- genomes[,1] # convert to vector
genomes

#i <- "Leptographium_lundbergii_S_CBS138716_GCA_001455505"

suppressMessages(
  for (i in genomes) {
    mcs <- NULL
    df <- read.table(paste0("kegg_output/",i,"/",i,".count_KEGG.txt"), sep="\t", header=1)
    kos <- df[,1]
    for (mid in mf) {
      mc <- ggkegg::module_completeness(ggkegg::module(mid, directory="mf"), query = kos)
      mcs <- c(mcs, mc$complete |> mean()) ## Mean of blocks
    }
    annos[[as.character(i)]] <- mcs
  }
)

warnings()


# We will next visualize the results using ComplexHeatmap and simplifyEnrichment. We will plot the word cloud of module description alongside heatmap by simplifyEnrichment, for determined clusters.

library(ComplexHeatmap)

## Make data.frame
hdf <- data.frame(annos, check.names=FALSE)
row.names(hdf) <- mf
hdf[is.na(hdf)] <- 0
hdf <- hdf[apply(hdf, 1, sum)!=0,] # remove modules with 0 completion for every genome
hdf <- as.matrix(hdf)

## Prepare for word cloud annotation
moddesc <- data.table::fread("https://rest.kegg.jp/list/module", header=FALSE)

## Obtain K-means clustering
km = kmeans(hdf, centers = 20)$cluster

gene_list <- split(row.names(hdf), km)
gene_list <- lapply(gene_list, function(x) {
  x[!is.na(x)]
})

annotList <- list()
for (i in names(gene_list)) {
  maps <- (moddesc |> dplyr::filter(V1 %in% gene_list[[i]]))$V2
  annotList[[i]] <-  maps
}


# generate heatmap
library(simplifyEnrichment)

#col_fun = circlize::colorRamp2(c(0, 0.5, 1), c(scales::muted("blue"), "white", scales::muted("red")))

col_fun = circlize::colorRamp2(c(0, 1),
                               c("white", scales::muted("blue")))

# rename heatmap columns to something simpler
column_labels <- genome
names(column_labels) = genomes
column_labels

ht1 <- Heatmap(hdf,
               name="Module\ncompleteness",
               show_column_names = TRUE,
               show_row_names = TRUE,
               row_names_side = "left",
               col=col_fun,
               row_split=km,
               cluster_columns = FALSE,
               #column_order = order(genomes),
               column_labels = column_labels[colnames(hdf)],
               column_names_gp = grid::gpar(fontsize = 7),
               row_names_gp = grid::gpar(fontsize = 6),
               row_title = NULL,
               width = unit(9, "cm"),
               height = unit(35, "cm"),
               heatmap_legend_param = list(
                 legend_direction = "vertical", 
                 legend_width = unit(5, "cm")
               ),
               rect_gp = gpar(col = "white", lwd = 0),
               border=TRUE,
               #column_names_max_height =unit(5,"cm")
              ) +
  rowAnnotation(
    keywords = simplifyEnrichment::anno_word_cloud(align_to = km,
                                                   term=annotList,
                                                   exclude_words=c("pathway",
                                                                   "degradation",
                                                                   "biosynthesis"),
                                                   max_words = 40,
                                                   fontsize_range = c(5,20))
  )

ht1

#htShiny() #interactive heatmap


pdf(file="~/Documents/2023_02_Sporothrix/results/figures/module_completeness.pdf", width=10, height=18)
draw(ht1, heatmap_legend_side = "left")
dev.off()

