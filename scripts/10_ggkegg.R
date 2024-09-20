#BiocManager::install("ggkegg")
library(ggkegg)
#devtools::install_github("noriakis/ggkegg")
library(tidyverse)
library(tidygraph)
library(ggfx)
library(ggraph)
library(igraph)
library(clusterProfiler)
library(kableExtra)
library(ComplexHeatmap)

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
## Prepare for word cloud annotation
moddesc <- data.table::fread("https://rest.kegg.jp/list/module", header=FALSE)
colnames(moddesc) = c("Module", "Description")
modrank <- data.table::fread("module_ranks.txt", header=TRUE) # how the modules are organized at https://www.genome.jp/brite/ko00002
moddesrank <- dplyr::left_join(modrank, moddesc, by = "Module")

## Make data.frame
hdf <- data.frame(annos, check.names=FALSE)
row.names(hdf) <- mf
hdf[is.na(hdf)] <- 0

#sort modules by rank
hdf.temp <- hdf %>% rownames_to_column(var = "Module")
hdf.sorted <- moddesrank %>% left_join(hdf.temp) %>% select(-cluster, -Rank_1, -Rank_2, -Rank_3, -Description) %>% column_to_rownames("Module")
hdf.sorted <- as.matrix(hdf.sorted)

#hdf <- hdf[apply(hdf, 1, sum)!=0,] # remove modules with 0 completion for every genome
hdf <- as.matrix(hdf)

## Obtain K-means clustering
km = kmeans(hdf, centers = 20)$cluster

gene_list <- split(row.names(hdf), km)
gene_list <- lapply(gene_list, function(x) {
  x[!is.na(x)]
})

annotList <- list()
for (i in names(gene_list)) {
  maps <- (moddesc |> dplyr::filter(Module %in% gene_list[[i]]))$Description
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

# manual clustering based on module category
cluster <- moddesrank$cluster
names(cluster) = moddesrank$Module

description_matrix <- moddesrank %>%
  column_to_rownames("Module") %>%
  select(Rank_3, Rank_2, Rank_1, Description) %>%
  as.matrix()

description_annotation <- rowAnnotation(
  "Module_category" = anno_text(description_matrix[,1], gp = gpar(fontsize = 6)),
  "Module_description" = anno_text(description_matrix[,4], gp = gpar(fontsize = 6)))

str(description_matrix)

ht2 <- Heatmap(hdf.sorted,
               name="Module\ncompleteness",
               show_column_names = TRUE,
               show_row_names = TRUE,
               row_names_side = "left",
               col=col_fun,
               row_split=cluster,
               cluster_columns = FALSE,
               #column_order = order(genomes),
               column_labels = column_labels[colnames(hdf.sorted)],
               column_names_gp = grid::gpar(fontsize = 7),
               cluster_rows = FALSE,
               #row_order = order(moddesrank$Module),
               #row_order = order(description_matrix),
               row_names_gp = grid::gpar(fontsize = 6),
               row_title = NULL,
               width = unit(9, "cm"),
               height = unit(100, "cm"),
               heatmap_legend_param = list(
                 legend_direction = "vertical", 
                 legend_width = unit(5, "cm")
               ),
               rect_gp = gpar(col = "white", lwd = 0),
               border=TRUE,
               right_annotation = description_annotation
)
  #rowAnnotation(
  #  keywords = simplifyEnrichment::anno_word_cloud(align_to = cluster,
  #                                                 term=annotList,
  #                                                 exclude_words=c("pathway",
  #                                                                 "degradation",
  #                                                                 "biosynthesis"),
  #                                                 max_words = 40,
  #                                                 fontsize_range = c(5,20))
 # )


ht2
pdf(file="~/Documents/2023_02_Sporothrix/results/figures/module_completeness_ranked.pdf", width=13, height=42)
draw(ht2, heatmap_legend_side = "left")
dev.off()






###########explore M000307########
module("M00307") |>
  module_text() |> ## return data.frame
  plot_module_text() ## wrapper function

mod <- module("M00307")
query <- c("K00161", "K00162", "K00382") #modules found in the genome
mod |>
  module_completeness(query) |>
  kableExtra::kable()



i <- "Leptographium_lundbergii_S_CBS138716_GCA_001455505"
i <- "Sporothrix_eucalyptigena_S_CBS139899_GCA_XXXXXXXXX"
mid <- "M00307"
annos <- list()
mcs <- NULL
df <- read.table(paste0("kegg_output/",i,"/",i,".count_KEGG.txt"), sep="\t", header=1)
kos <- df[,1]
mc <- ggkegg::module_completeness(ggkegg::module(mid, directory="mf"), query = kos) #the kos of that genome in that module
mcs <- c(mcs, mc$complete |> mean()) ## Mean of blocks
annos[[as.character(i)]] <- mcs
annos


query <- c("K00163","K00161","K00162","K00627","K00382","K13997","K00169","K00170","K00171","K00172","K00189","K03737")

#mod <- mc$block

#query <- gsub("\\+","&",gsub(" ", "&", gsub(",", "|", mod)))
#query <- gsub("\\-", "+0*", mod)
#query

df %>% filter(df[,1] %in% c(query)) #which KOs are present in the given genome and module

mod



###########explore M00897########
module("M00897") |>
  module_text("2") |> ## return data.frame
  plot_module_text() ## wrapper function

mod <- module("M00897")
module_kos <- c("K03146", "K03147", "K14153", "K22911", "K00949") # kos in the module

i <- "Sporothrix_eucalyptigena_S_CBS139899_GCA_XXXXXXXXX"
mid <- "M00897"
annos <- list()
mcs <- NULL
df <- read.table(paste0("kegg_output/",i,"/",i,".count_KEGG.txt"), sep="\t", header=1)
kos <- df[,1]
mc <- ggkegg::module_completeness(ggkegg::module(mid, directory="mf"), query = kos, name="2") #the kos of that genome in that module
mcs <- c(mcs, mc$complete |> mean()) ## Mean of blocks
annos[[as.character(i)]] <- mcs
annos


#mod <- mc$block

#query <- gsub("\\+","&",gsub(" ", "&", gsub(",", "|", mod)))
#query <- gsub("\\-", "+0*", mod)
#query

query <- df %>% filter(.[[1]] %in% c(module_kos)) %>% pull(paste0(i)) #which KOs are present in the given genome and module

mod |>
  module_completeness(query) |>
  kableExtra::kable()


###########explore M00009########
module("M00009") |>
  module_text() |> ## return data.frame
  plot_module_text() ## wrapper function

mod <- module("M00009")
mod@definition



module_kos <- c("K01647", "K05942", "K01681", "K01682", "K00031", "K00030", "K00164", "K00658", "K01616", "K00382", "K00174", "K00175", "K00177", "K00176", "K01902", "K01903", "K01899", "K01900", "K18118", "K00234", "K00235", "K00236", "K00237", "K25801", "K00239", "K00240", "K00241", "K00242", "K18859", "K18860", "K00244", "K00245", "K00246", "K00247", "K01676", "K01679", "K01677", "K01678", "K00026", "K00025", "K00024", "K00116") # kos in the module

i <- "Sporothrix_eucalyptigena_S_CBS139899_GCA_XXXXXXXXX"
mid <- "M00009"
annos <- list()
mcs <- NULL
df <- read.table(paste0("kegg_output/",i,"/",i,".count_KEGG.txt"), sep="\t", header=1)
kos <- df[,1]
mc <- ggkegg::module_completeness(ggkegg::module(mid, directory="mf"), query = kos) #the kos of that genome in that module
mcs <- c(mcs, mc$complete |> mean()) ## Mean of blocks
annos[[as.character(i)]] <- mcs
annos


#mod <- mc$block

#query <- gsub("\\+","&",gsub(" ", "&", gsub(",", "|", mod)))
#query <- gsub("\\-", "+0*", mod)
#query

query <- df %>% filter(.[[1]] %in% c(module_kos)) %>% pull(paste0(i)) #which KOs are present in the given genome and module

mod |>
  module_completeness(query) |>
  kableExtra::kable()

###########explore M00897########
module("M00897") |>
  module_text() |> ## return data.frame
  plot_module_text() ## wrapper function

mod <- module("M00897")
module_kos <- c("K03146", "K03147", "K14153", "K22911", "K00949") # kos in the module

i <- "Sporothrix_eucalyptigena_S_CBS139899_GCA_XXXXXXXXX"
mid <- "M00897"
annos <- list()
mcs <- NULL
df <- read.table(paste0("kegg_output/",i,"/",i,".count_KEGG.txt"), sep="\t", header=1)
kos <- df[,1]
mc <- ggkegg::module_completeness(ggkegg::module(mid, directory="mf"), query = kos) #the kos of that genome in that module
mcs <- c(mcs, mc$complete |> mean()) ## Mean of blocks
annos[[as.character(i)]] <- mcs
annos


#mod <- mc$block

#query <- gsub("\\+","&",gsub(" ", "&", gsub(",", "|", mod)))
#query <- gsub("\\-", "+0*", mod)
#query

query <- df %>% filter(.[[1]] %in% c(module_kos)) %>% pull(paste0(i)) #which KOs are present in the given genome and module

mod |>
  module_completeness(query) |>
  kableExtra::kable()


###########explore M00899 thiamine salvage#######
module("M00899") |>
  module_text() |> ## return data.frame
  plot_module_text() ## wrapper function

mod <- module("M00899")
mod@definition_raw


module_kos <- c("K01647", "K05942", "K01681", "K01682", "K00031", "K00030", "K00164", "K00658", "K01616", "K00382", "K00174", "K00175", "K00177", "K00176", "K01902", "K01903", "K01899", "K01900", "K18118", "K00234", "K00235", "K00236", "K00237", "K25801", "K00239", "K00240", "K00241", "K00242", "K18859", "K18860", "K00244", "K00245", "K00246", "K00247", "K01676", "K01679", "K01677", "K01678", "K00026", "K00025", "K00024", "K00116") # kos in the module

i <- "Sporothrix_eucalyptigena_S_CBS139899_GCA_XXXXXXXXX"
mid <- "M00009"
annos <- list()
mcs <- NULL
df <- read.table(paste0("kegg_output/",i,"/",i,".count_KEGG.txt"), sep="\t", header=1)
kos <- df[,1]
mc <- ggkegg::module_completeness(ggkegg::module(mid, directory="mf"), query = kos) #the kos of that genome in that module
mcs <- c(mcs, mc$complete |> mean()) ## Mean of blocks
annos[[as.character(i)]] <- mcs
annos


#mod <- mc$block

#query <- gsub("\\+","&",gsub(" ", "&", gsub(",", "|", mod)))
#query <- gsub("\\-", "+0*", mod)
#query

query <- df %>% filter(.[[1]] %in% c(module_kos)) %>% pull(paste0(i)) #which KOs are present in the given genome and module

mod <- module("M00009")
mod |>
  module_completeness(query) |>
  kableExtra::kable()



module("M00577") |>
  module_text() |> ## return data.frame
  plot_module_text() ## wrapper function


mod <- module("M00899")
mod <- module("M00009")

mod@definitions
mod@definitions[["1"]]
mod@definitions[["2"]]

## Extract definition 2
mod |>
  module_text("2") |> ## return data.frame
  plot_module_text() ## wrapper function
  