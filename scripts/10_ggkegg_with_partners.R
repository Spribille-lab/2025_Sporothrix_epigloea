#BiocManager::install("ggkegg")
library(ggkegg)
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

genome <- c("Sporothrix_euskadiensis_VPRI43754",
            "Sporothrix_pseudoabietina_VPRI43531",
            "Sporothrix_variecibatus_CBS121961",
            "Sporothrix_protearum_CBS116654",
            "Sporothrix_humicola_CBS118129",
            "Sporothrix_pallida_CBS13156",
            "Sporothrix_mexicana_CBS120341",
            "Sporothrix_brasiliensis_5110",
            "Sporothrix_schenckii_1099",
            "Sporothrix_globosa_CBS120340",
            "Sporothrix_luriei_CBS93772",
            "Sporothrix_phasma_CBS119721",
            "Sporothrix_dimorphospora_CBS55374",
            "Sporothrix_inflata_CBS23968",
            "Sporothrix_bragantina_CBS47491",
            "Sporothrix_curviconia_CBS95973",
            "Sporothrix_thermara_CBS139747MAG",
            "Sporothrix_epigloea_CBS57363",
            "Sporothrix_epigloea_TF4163",
            "Sporothrix_epigloea_CBS119000",
            "Sporothrix_eucalyptigena_CBS139899",
            "Sporothrix_eucalyptigena_CBS140593",
            "Sporothrix_nigrograna_VPRI43755",
            "Ophiostoma_novoulmi_H327",
            "Ophiostoma_ips_VPRI43529",
            "Sporothrix_brunneoviolacea_CBS124561",
            "Ophiostoma_fasciatum_VPRI43845",
            "Sporothrix_insectorum_RCEF264",
            "Leptographium_lundbergii_CBS138716",
            "Tremella_yokohamensis_CBS18117",
            "Annulohypoxylon_annulatum_CBS149473")
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
genomes <- read.delim("../genomes_with_partners.txt")
genomes <- genomes[,1] # convert to vector
genomes

# module annotations
moddesc <- data.table::fread("https://rest.kegg.jp/list/module", header=FALSE)
colnames(moddesc) = c("Module", "Description")
modrank <- data.table::fread("module_ranks.txt", header=TRUE) # how the modules are organized at https://www.genome.jp/brite/ko00002
moddesrank <- dplyr::left_join(modrank, moddesc, by = "Module")

#i <- "Leptographium_lundbergii_S_CBS138716_GCA_001455505"


#########################
#function

#genomes=c("Sporothrix_euskadiensis_S_VPRI43754_GCA_019925375","Sporothrix_pseudoabietina_S_VPRI43531_GCA_019925295") #test with only two 

module_completeness_by_definition <- function(definition) {
  suppressMessages(
    for (i in genomes) {
      mcs <- NULL
      df <- read.table(paste0("kegg_output/",i,"/",i,".count_KEGG.txt"), sep="\t", header=1)
      kos <- df[,1]
      for (mid in mf) {
        mc <- ggkegg::module_completeness(ggkegg::module(mid, directory="mf"), query = kos, name = definition)
        mcs <- c(mcs, mc$complete |> mean()) ## Mean of blocks
      }
      annos[[as.character(i)]] <- mcs
    }
  )
  
  # Make data.frame
  hdf <- data.frame(annos, check.names=FALSE)
  row.names(hdf) <- mf
  hdf[is.na(hdf)] <- 0
  
  #sort modules by rank
  hdf.temp <- hdf %>% rownames_to_column(var = "Module")
  hdf.sorted <- moddesrank %>% left_join(hdf.temp) %>% select(-cluster, -Rank_1, -Rank_2, -Rank_3, -Description)
  hdf.sorted$definition <- definition #create a column with the definition name
  hdf.sorted <- hdf.sorted %>% relocate(definition) %>% unite("Module_def", c("Module", "definition"), sep = "_", remove = FALSE)
  hdf.sorted$order <- c(1:nrow(hdf.sorted))
  hdf.sorted <- hdf.sorted %>% relocate(order)
  return(hdf.sorted)
}

definitions=c("1","2","3","4") # there are a maximum of 4 definitions per module as far as I could tell

l<-lapply(definitions, module_completeness_by_definition) #apply the function to each definition
hdf_combined <- do.call(rbind, l)
str(hdf_combined)


#join descriptions with hdf

hdf_combined_desc <- left_join(hdf_combined, moddesrank, by = join_by(Module))
hdf_combined_desc <- hdf_combined_desc %>% dplyr::arrange(order)

# set up data for heatmap in a matrix
hdf_combined <- hdf_combined_desc %>% dplyr::select(-order, -definition, -Module, -cluster, -Rank_1, -Rank_2, -Rank_3, -Description)
hdf_combined <- hdf_combined %>% column_to_rownames("Module_def")
hdf_combined <- hdf_combined[apply(hdf_combined, 1, sum)!=0,]

# set up descriptions for heatmap
annotations <- hdf_combined_desc %>% dplyr::select(Module_def, cluster, Rank_1, Rank_2, Rank_3, Description)
annotations <- annotations %>% column_to_rownames("Module_def")
annotations <- annotations %>% filter(row.names(annotations) %in% row.names(hdf_combined)) #remove description for modules that aren't included in the heatmap



# generate heatmap
col_fun = circlize::colorRamp2(c(0, 1),
                               c("white", scales::muted("blue")))

# rename heatmap columns to something simpler
column_labels <- genome
names(column_labels) = genomes
column_labels

# manual clustering based on module category
cluster <- annotations$cluster
names(cluster) = annotations$Module_def
cluster

cluster_col <- c(rep(1,17), rep(2,3), rep(3,9), rep(4,2))
names(cluster_col) = genomes
cluster_col

description_matrix <- annotations %>%
  select(Rank_3, Rank_2, Rank_1, Description) %>%
  as.matrix()

description_annotation <- rowAnnotation(
  "Module_category" = anno_text(description_matrix[,1], gp = gpar(fontsize = 8)),
  "Module_description" = anno_text(description_matrix[,4], gp = gpar(fontsize = 8)))

str(description_matrix)

hdf_combined <- as.matrix(hdf_combined)

ht2 <- Heatmap(hdf_combined,
               name="Module\ncompleteness",
               show_column_names = TRUE,
               show_row_names = TRUE,
               row_names_side = "left",
               col=col_fun,
               row_split=cluster,
               column_split=cluster_col,
               cluster_columns = FALSE,
               #column_order = order(genomes),
               column_labels = column_labels[colnames(hdf_combined)],
               column_names_gp = grid::gpar(fontsize = 7),
               cluster_rows = FALSE,
               #row_order = order(moddesrank$Module),
               #row_order = order(description_matrix),
               row_names_gp = grid::gpar(fontsize = 8),
               row_title = NULL,
               column_title = NULL,
               width = unit(9, "cm"),
               height = unit(60, "cm"),
               heatmap_legend_param = list(
                 legend_direction = "vertical", 
                 legend_width = unit(5, "cm")
               ),
               rect_gp = gpar(col = "white", lwd = 0),
               border=TRUE,
               right_annotation = description_annotation
)
ht2
pdf(file="~/Documents/2023_02_Sporothrix/results/figures/module_completeness_with_partners.pdf",
    width=13, height=30)
draw(ht2, heatmap_legend_side = "left")
dev.off()
