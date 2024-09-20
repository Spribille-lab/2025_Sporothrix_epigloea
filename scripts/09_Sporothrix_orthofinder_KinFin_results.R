library(tidyverse)
library("UpSetR")
# library(seqRFLP)
library(data.table)
library(ComplexHeatmap)
library(phylotools)
library(readxl)

sessionInfo()

setwd("~/Documents/2023_02_Sporothrix/09_orthology/OrthoFinder/Results_Aug07/")

# read orthofinder output

#file_list <- list.files(path="../../../orthofinder_input/")
genome_list = c("Sporothrix_bragantina_CBS47491",
                "Sporothrix_brasiliensis_5110",
                "Sporothrix_curviconia_CBS95973",
                "Sporothrix_dimorphospora_CBS55374",
                "Sporothrix_epigloea_CBS119000",
                "Sporothrix_epigloea_CBS57363",
                "Sporothrix_epigloea_TF4163",
                "Sporothrix_eucalyptigena_CBS139899",
                "Sporothrix_eucalyptigena_CBS140593",
                "Sporothrix_euskadiensis_VPRI43754",
                "Sporothrix_globosa_CBS120340",
                "Sporothrix_humicola_CBS118129",
                "Sporothrix_inflata_CBS23968",
                "Sporothrix_luriei_CBS93772",
                "Sporothrix_mexicana_CBS120341",
                "Sporothrix_nigrograna_VPRI43755",
                "Sporothrix_pallida_CBS13156",
                "Sporothrix_phasma_CBS119721",
                "Sporothrix_protearum_CBS116654",
                "Sporothrix_pseudoabietina_VPRI43531",
                "Sporothrix_schenckii_1099",
                "Sporothrix_thermara_CBS139747MAG",
                "Sporothrix_variecibatus_CBS121961")

# grab the Orthogroups.tsv file from orthofinder, write a table for each genome, write each vector to file.  I can't quite remember why I wrote everythng to a file and then imported it again.
for (i in 1:length(genome_list)){
  orthogroups <- read.delim("Orthogroups/Orthogroups.tsv", stringsAsFactors = F) %>%
  select(Orthogroup, paste(genome_list[i]))
colnames(orthogroups) = c("Orthogroup", "Protein_ID")
orthogroups <- orthogroups %>%
  separate_rows(Protein_ID, sep = ", ") %>%
  filter(Protein_ID != "")
write.table(orthogroups, paste("Custom_results/",genome_list[i],".orthogroups.txt", sep = ""), sep='\t', quote = F, row.names = F)
orthogroup_v <- pull(orthogroups, Orthogroup)
saveRDS(orthogroup_v, file = paste("Custom_results/",genome_list[i],".orthogroups.rds", sep=""))
orthogroup_unique_v = unique(orthogroup_v)
saveRDS(orthogroup_unique_v, file = paste("Custom_results/",genome_list[i],".unique.orthogroups.rds", sep=""))
}

# make an object for each vector of orthogroups. I wrote this before I knew how to write functions.  Don't judge me.
Sporothrix_bragantina_CBS47491 <- readRDS("Custom_results/Sporothrix_bragantina_CBS47491.unique.orthogroups.rds")
Sporothrix_brasiliensis_5110 <- readRDS("Custom_results/Sporothrix_brasiliensis_5110.unique.orthogroups.rds")
Sporothrix_curviconia_CBS95973 <- readRDS("Custom_results/Sporothrix_curviconia_CBS95973.unique.orthogroups.rds")
Sporothrix_dimorphospora_CBS55374 <- readRDS("Custom_results/Sporothrix_dimorphospora_CBS55374.unique.orthogroups.rds")
Sporothrix_epigloea_CBS119000 <- readRDS("Custom_results/Sporothrix_epigloea_CBS119000.unique.orthogroups.rds")
Sporothrix_epigloea_TF4163 <- readRDS("Custom_results/Sporothrix_epigloea_TF4163.unique.orthogroups.rds")
Sporothrix_epigloea_CBS57363 <- readRDS("Custom_results/Sporothrix_epigloea_CBS57363.unique.orthogroups.rds")
Sporothrix_eucalyptigena_CBS139899 <- readRDS("Custom_results/Sporothrix_eucalyptigena_CBS139899.unique.orthogroups.rds")
Sporothrix_eucalyptigena_CBS140593 <- readRDS("Custom_results/Sporothrix_eucalyptigena_CBS140593.unique.orthogroups.rds")
Sporothrix_euskadiensis_VPRI43754 <- readRDS("Custom_results/Sporothrix_euskadiensis_VPRI43754.unique.orthogroups.rds")
Sporothrix_globosa_CBS120340 <- readRDS("Custom_results/Sporothrix_globosa_CBS120340.unique.orthogroups.rds")
Sporothrix_humicola_CBS118129 <- readRDS("Custom_results/Sporothrix_humicola_CBS118129.unique.orthogroups.rds")
Sporothrix_inflata_CBS23968 <- readRDS("Custom_results/Sporothrix_inflata_CBS23968.unique.orthogroups.rds")
Sporothrix_luriei_CBS93772 <- readRDS("Custom_results/Sporothrix_luriei_CBS93772.unique.orthogroups.rds")
Sporothrix_mexicana_CBS120341 <- readRDS("Custom_results/Sporothrix_mexicana_CBS120341.unique.orthogroups.rds")
Sporothrix_nigrograna_VPRI43755 <- readRDS("Custom_results/Sporothrix_nigrograna_VPRI43755.unique.orthogroups.rds")
Sporothrix_pallida_CBS13156 <- readRDS("Custom_results/Sporothrix_pallida_CBS13156.unique.orthogroups.rds")
Sporothrix_phasma_CBS119721 <- readRDS("Custom_results/Sporothrix_phasma_CBS119721.unique.orthogroups.rds")
Sporothrix_protearum_CBS116654 <- readRDS("Custom_results/Sporothrix_protearum_CBS116654.unique.orthogroups.rds")
Sporothrix_pseudoabietina_VPRI43531 <- readRDS("Custom_results/Sporothrix_pseudoabietina_VPRI43531.unique.orthogroups.rds")
Sporothrix_schenckii_1099 <- readRDS("Custom_results/Sporothrix_schenckii_1099.unique.orthogroups.rds")
Sporothrix_thermara_CBS139747MAG <- readRDS("Custom_results/Sporothrix_thermara_CBS139747MAG.unique.orthogroups.rds") 
Sporothrix_variecibatus_CBS121961 <- readRDS("Custom_results/Sporothrix_variecibatus_CBS121961.unique.orthogroups.rds") 

##################################

# all sporothrix genomes

sporothrix_list <- list(
                      S_epigloea_TF4163 = Sporothrix_epigloea_TF4163,
                      S_epigloea_CBS119000 = Sporothrix_epigloea_CBS119000,
                      S_epigloea_CBS57363 = Sporothrix_epigloea_CBS57363,
                      S_curviconia = Sporothrix_curviconia_CBS95973,
                      S_nigrograna = Sporothrix_nigrograna_VPRI43755,
                      S_brasiliensis = Sporothrix_brasiliensis_5110,
                      S_dimorphospora = Sporothrix_dimorphospora_CBS55374,
                      S_euskadiensis = Sporothrix_euskadiensis_VPRI43754,
                      S_globosa = Sporothrix_globosa_CBS120340,
                      S_humicola = Sporothrix_humicola_CBS118129,
                      S_inflata =  Sporothrix_inflata_CBS23968,
                      S_luriei = Sporothrix_luriei_CBS93772,
                      S_mexicana = Sporothrix_mexicana_CBS120341,
                      S_pallida = Sporothrix_pallida_CBS13156,
                      S_phasma = Sporothrix_phasma_CBS119721,
                      S_protearum = Sporothrix_protearum_CBS116654,
                      S_pseudoabietina = Sporothrix_pseudoabietina_VPRI43531,
                      S_schenckii = Sporothrix_schenckii_1099,
                      S_variecibatus = Sporothrix_variecibatus_CBS121961,
                      S_thermara = Sporothrix_thermara_CBS139747MAG,
                      S_eucalyptigena_CBS139899 = Sporothrix_eucalyptigena_CBS139899,
                      S_eucalyptigena_CBS140593 = Sporothrix_eucalyptigena_CBS140593,
                      S_bragantina = Sporothrix_bragantina_CBS47491)

# make dataframe for metadata
sets = c("S_epigloea_TF4163","S_epigloea_CBS119000","S_epigloea_CBS57363","S_curviconia","S_nigrograna","S_brasiliensis","S_dimorphospora","S_euskadiensis","S_globosa", "S_humicola", "S_inflata", "S_luriei","S_mexicana","S_pallida","S_phasma","S_protearum","S_pseudoabietina","S_schenckii","S_variecibatus", "S_thermara", "S_eucalyptigena_CBS139899", "S_eucalyptigena_CBS140593", "S_bragantina")
group = c("epigloea", "epigloea", "epigloea", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other")
metadata <- as.data.frame(cbind(sets, group))
  
  
pdf(file="~/Documents/2023_02_Sporothrix/09_orthology/OrthoFinder/Results_Aug07/figures/orthofinder_upset_core_sporothrix.pdf",
    width=15, height=9)
help("upset")

upset(fromList(sporothrix_list),
      nsets = 23,
      keep.order = TRUE,
      order.by = "freq",
      mainbar.y.label = "Intersection size",
      sets.x.label = "Unique orthogroups",
      scale.intersections = "identity",
      #mainbar.y.max = 500,
      set.metadata = list(data = metadata,
                          plots = list(list(type = "matrix_rows",
                                            column = "group",
                                            colors = c(epigloea = "azure3", other = "white"),
                                            alpha = 0.5))),
      #intersections = list(list("S_insectorum","S_brunneoviolacea","S_nigrograna","S_epigloea_CBS119000","S_epigloea_CBS57363","S_epigloea_TF4163","S_inflata","S_dimorphospora","S_phasma","S_luriei","S_globosa","S_schenckii","S_brasiliensis","S_mexicana","S_pallida","S_humicola","S_protearum","S_variecibatus","S_euskadiensis","S_pseudoabietina"),
       #                    list("S_insectorum","S_brunneoviolacea","S_nigrograna","S_inflata","S_dimorphospora","S_phasma","S_luriei","S_globosa","S_schenckii","S_brasiliensis","S_mexicana","S_pallida","S_humicola","S_protearum","S_variecibatus","S_euskadiensis","S_pseudoabietina"),
      #                     list("S_brunneoviolacea","S_nigrograna","S_epigloea_CBS119000","S_epigloea_CBS57363","S_epigloea_TF4163","S_inflata","S_dimorphospora","S_phasma","S_luriei","S_globosa","S_schenckii","S_brasiliensis","S_mexicana","S_pallida","S_humicola","S_protearum","S_variecibatus","S_euskadiensis","S_pseudoabietina"),
      #                     list("S_pallida","S_humicola"),
      #                     list("S_insectorum","S_brunneoviolacea","S_nigrograna","S_inflata","S_dimorphospora","S_luriei","S_globosa","S_schenckii","S_brasiliensis","S_mexicana","S_pallida","S_humicola","S_protearum","S_variecibatus","S_euskadiensis","S_pseudoabietina"),
      #                     list("S_mexicana","S_pallida","S_humicola"),
      #                     list("S_insectorum","S_brunneoviolacea"),
      #                     list("S_brunneoviolacea","S_nigrograna","S_inflata","S_dimorphospora","S_phasma","S_luriei","S_globosa","S_schenckii","S_brasiliensis","S_mexicana","S_pallida","S_humicola","S_protearum","S_variecibatus","S_euskadiensis","S_pseudoabietina"),
      #                     list("S_insectorum","S_brunneoviolacea","S_epigloea_CBS119000","S_epigloea_CBS57363","S_epigloea_TF4163","S_inflata","S_dimorphospora","S_phasma","S_luriei","S_globosa","S_schenckii","S_brasiliensis","S_mexicana","S_pallida","S_humicola","S_protearum","S_variecibatus","S_euskadiensis","S_pseudoabietina"),
       #                    list("S_insectorum","S_brunneoviolacea","S_nigrograna","S_epigloea_CBS119000","S_epigloea_CBS57363","S_epigloea_TF4163","S_inflata","S_dimorphospora","S_phasma","S_luriei","S_globosa","S_schenckii","S_brasiliensis","S_mexicana","S_pallida","S_humicola","S_variecibatus","S_euskadiensis","S_pseudoabietina"),
      #                     list("S_insectorum","S_brunneoviolacea","S_nigrograna","S_epigloea_CBS119000","S_epigloea_CBS57363","S_epigloea_TF4163","S_inflata","S_dimorphospora","S_phasma","S_globosa","S_schenckii","S_brasiliensis","S_mexicana","S_pallida","S_humicola","S_protearum","S_variecibatus","S_euskadiensis","S_pseudoabietina"),
       #                    list("S_insectorum","S_brunneoviolacea","S_inflata","S_dimorphospora","S_luriei","S_globosa","S_schenckii","S_brasiliensis","S_mexicana","S_pallida","S_humicola","S_protearum","S_variecibatus"),
      #                     list("S_insectorum","S_nigrograna","S_epigloea_CBS119000","S_epigloea_CBS57363","S_epigloea_TF4163","S_inflata","S_dimorphospora","S_phasma","S_luriei","S_globosa","S_schenckii","S_brasiliensis","S_mexicana","S_pallida","S_humicola","S_protearum","S_variecibatus","S_euskadiensis","S_pseudoabietina"),
       #                    list("S_insectorum","S_brunneoviolacea","S_inflata","S_dimorphospora","S_phasma","S_luriei","S_globosa","S_schenckii","S_brasiliensis","S_mexicana","S_pallida","S_humicola","S_protearum","S_variecibatus","S_euskadiensis","S_pseudoabietina"),
        #                   list("S_insectorum","S_brunneoviolacea","S_inflata","S_dimorphospora","S_luriei","S_globosa","S_schenckii","S_brasiliensis","S_mexicana","S_pallida","S_humicola","S_variecibatus"),
       #                   list("S_insectorum","S_brunneoviolacea","S_nigrograna","S_epigloea_CBS57363","S_epigloea_TF4163","S_inflata","S_dimorphospora","S_phasma","S_luriei","S_globosa","S_schenckii","S_brasiliensis","S_mexicana","S_pallida","S_humicola","S_protearum","S_variecibatus","S_euskadiensis","S_pseudoabietina"),
       #                    list("S_insectorum","S_brunneoviolacea","S_inflata","S_dimorphospora","S_luriei","S_globosa","S_schenckii","S_brasiliensis","S_mexicana","S_pallida","S_humicola","S_protearum","S_variecibatus","S_euskadiensis","S_pseudoabietina"),
       #                    list("S_epigloea_CBS57363","S_epigloea_TF4163"),
       #                    list("S_insectorum","S_brunneoviolacea","S_nigrograna","S_epigloea_CBS119000","S_epigloea_CBS57363","S_epigloea_TF4163","S_inflata","S_dimorphospora","S_luriei","S_globosa","S_schenckii","S_brasiliensis","S_mexicana","S_pallida","S_humicola","S_protearum","S_variecibatus","S_euskadiensis","S_pseudoabietina"),
       #                    list("S_brunneoviolacea","S_inflata","S_dimorphospora"),
        #                   list("S_brunneoviolacea","S_nigrograna","S_inflata","S_dimorphospora","S_luriei","S_globosa","S_schenckii","S_brasiliensis","S_mexicana","S_pallida","S_humicola","S_protearum","S_variecibatus","S_euskadiensis","S_pseudoabietina"),
        #                   list("S_inflata","S_dimorphospora"),
         #                  list("S_epigloea_CBS119000","S_epigloea_CBS57363","S_epigloea_TF4163"),
         #                  list("S_brunneoviolacea","S_dimorphospora"),
        #                   list("S_brunneoviolacea","S_inflata"),
        #                   list("S_brunneoviolacea","S_inflata","S_dimorphospora","S_phasma","S_luriei","S_globosa","S_schenckii","S_brasiliensis","S_mexicana","S_pallida","S_humicola","S_protearum","S_variecibatus","S_euskadiensis","S_pseudoabietina"),
        #                   list("S_insectorum","S_brunneoviolacea","S_inflata","S_dimorphospora","S_luriei","S_globosa","S_schenckii","S_brasiliensis","S_mexicana","S_pallida","S_humicola","S_variecibatus","S_euskadiensis","S_pseudoabietina"),
        #                   list("S_insectorum","S_brunneoviolacea","S_nigrograna","S_inflata","S_dimorphospora","S_luriei","S_globosa","S_schenckii","S_brasiliensis","S_mexicana","S_pallida","S_humicola","S_variecibatus","S_euskadiensis","S_pseudoabietina"),
        #                   list("S_epigloea_CBS119000"),
        #                   list("S_brunneoviolacea","S_inflata","S_dimorphospora","S_luriei","S_globosa","S_schenckii","S_brasiliensis","S_mexicana","S_pallida","S_humicola","S_protearum","S_variecibatus","S_euskadiensis","S_pseudoabietina"),
         #                  list("S_brunneoviolacea","S_inflata","S_dimorphospora","S_luriei","S_globosa","S_schenckii","S_brasiliensis","S_mexicana","S_pallida","S_humicola","S_variecibatus","S_euskadiensis","S_pseudoabietina"),
         #                  list("S_euskadiensis","S_pseudoabietina"),
         #                  list("S_brunneoviolacea","S_inflata","S_dimorphospora","S_mexicana","S_pallida","S_humicola","S_variecibatus"),
         #                  list("S_insectorum"),
         #                  list("S_insectorum","S_brunneoviolacea","S_inflata","S_dimorphospora","S_mexicana","S_pallida","S_humicola","S_variecibatus"),
         #                  list("S_variecibatus","S_euskadiensis","S_pseudoabietina"),
         #                  list("S_nigrograna","S_inflata","S_dimorphospora","S_phasma","S_luriei","S_globosa","S_schenckii","S_brasiliensis","S_mexicana","S_pallida","S_humicola","S_protearum","S_variecibatus","S_euskadiensis","S_pseudoabietina"),
         #                  list("S_nigrograna","S_epigloea_CBS119000","S_epigloea_CBS57363","S_epigloea_TF4163","S_inflata","S_dimorphospora","S_phasma","S_luriei","S_globosa","S_schenckii","S_brasiliensis","S_mexicana","S_pallida","S_humicola","S_protearum","S_variecibatus","S_euskadiensis","S_pseudoabietina"),
          #                 list("S_protearum","S_variecibatus","S_euskadiensis","S_pseudoabietina"),
         #                  list("S_brunneoviolacea")),
      sets = c("S_nigrograna",
               "S_eucalyptigena_CBS140593",
               "S_eucalyptigena_CBS139899",
               "S_curviconia",
               "S_bragantina",
               "S_thermara",
               "S_epigloea_CBS119000",
               "S_epigloea_TF4163",
               "S_epigloea_CBS57363",
               "S_inflata",
               "S_dimorphospora",
               "S_phasma",
               "S_luriei",
               "S_globosa",
               "S_schenckii",
               "S_brasiliensis",
               "S_mexicana",
               "S_humicola",
               "S_pallida",
               "S_protearum",
               "S_variecibatus",
               "S_euskadiensis",
               "S_pseudoabietina"
               ))

dev.off()

##### use ComplexHeatmap to extract combination sets -------

# convert from list to the binary matrix (ComplexHeatmaps)
sporothrix_mat <- list_to_matrix(sporothrix_list)

# generate the combination matrix and calculate the size of the sets and the combination sets (ComplexHeatmaps)

set.seed(123)
sporothrix_comb_mat = make_comb_mat(list_to_matrix(sporothrix_list))
sporothrix_comb_mat_t = t(sporothrix_comb_mat)
sporothrix_comb_mat
sporothrix_comb_mat[2:2377] #drop the largest combination set
str(sporothrix_comb_mat)

set_name(sporothrix_comb_mat)
comb_name(sporothrix_comb_mat) 
set_size(sporothrix_comb_mat) #The set sizes.
comb_size(sporothrix_comb_mat) #The combination set sizes

# extract orthogroups that are missing in epigloea and present in at least 19 other sporothrix genomes
OG1 <- extract_comb(sporothrix_comb_mat, "00011111111111111111111")
OG2 <- extract_comb(sporothrix_comb_mat, "00001111111111111111111")
OG3 <- extract_comb(sporothrix_comb_mat, "00010111111111111111111")
OG4 <- extract_comb(sporothrix_comb_mat, "00011011111111111111111")
OG5 <- extract_comb(sporothrix_comb_mat, "00011101111111111111111")
OG6 <- extract_comb(sporothrix_comb_mat, "00011110111111111111111")
OG7 <- extract_comb(sporothrix_comb_mat, "00011111011111111111111")
#OG8 <- extract_comb(sporothrix_comb_mat, "00011111101111111111111")
OG9 <- extract_comb(sporothrix_comb_mat, "00011111110111111111111")
OG10 <- extract_comb(sporothrix_comb_mat, "00011111111011111111111")
#OG11 <- extract_comb(sporothrix_comb_mat, "00011111111101111111111")
#OG12 <- extract_comb(sporothrix_comb_mat, "00011111111110111111111")
OG13 <- extract_comb(sporothrix_comb_mat, "00011111111111011111111")
OG14 <- extract_comb(sporothrix_comb_mat, "00011111111111101111111")
OG15 <- extract_comb(sporothrix_comb_mat, "00011111111111110111111")
OG16 <- extract_comb(sporothrix_comb_mat, "00011111111111111011111")
OG17 <- extract_comb(sporothrix_comb_mat, "00011111111111111101111")
OG18 <- extract_comb(sporothrix_comb_mat, "00011111111111111110111")
#OG19 <- extract_comb(sporothrix_comb_mat, "00011111111111111111011")
OG20 <- extract_comb(sporothrix_comb_mat, "00011111111111111111101")
OG21 <- extract_comb(sporothrix_comb_mat, "00011111111111111111110")


OG23 <- extract_comb(sporothrix_comb_mat, "11100000000000000000000")
OG24 <- extract_comb(sporothrix_comb_mat, "10100000000000000000000")
OG25 <- extract_comb(sporothrix_comb_mat, "01000000000000000000000")


#lost_by_epigloea <- c(OG1, OG2, OG3, OG4, OG5, OG6, OG7, OG9, OG10, OG13, OG14, OG15, OG16, OG17, OG18, OG20, OG21)
lost_by_epigloea <- OG1

lost_by_epigloea_df <- as.data.frame(lost_by_epigloea, stringsAsFactors = FALSE)
write.table(lost_by_epigloea_df, "gained_and_lost_strict/epigloea_orthogroups_lost.txt", row.names = FALSE, sep = "\t", col.names = FALSE, quote = FALSE)

#gained_by_epigloea <- c(OG23, OG24, OG25)
gained_by_epigloea <- OG23
gained_by_epigloea_df <- as.data.frame(gained_by_epigloea, stringsAsFactors = FALSE)
write.table(gained_by_epigloea_df, "gained_and_lost_strict/epigloea_orthogroups_gained.txt", row.names = FALSE, sep = "\t", col.names = FALSE, quote = FALSE)



########### make dataframe of orthgroups lost and gained --------
epigloea_list = c("Sporothrix_epigloea_S_CBS119000_GCA_943908295",
                "Sporothrix_epigloea_S_CBS57363_GCA_943900835",
                "Sporothrix_epigloea_S_TF4163MAG_GCA_944036445")

epigloea_annotations <- data.frame()
for (i in 1:length(epigloea_list)){
    temp_data <- read.delim(paste("~/Documents/2023_02_Sporothrix/10_annotations/annotation_tables_merged/",epigloea_list[i],".annotations.merged.txt", sep = ""), stringsAsFactors = F)
    temp_data <- temp_data %>%
      #select(Transcript_ID,Orthogroup,Translation) %>%
      add_column(genome = epigloea_list[i], .before = "Transcript_ID") %>%
      #dplyr::mutate(genome = str_replace(genome, "(\\S*)(.annotations.txt)", "\\1")) %>%
      dplyr::mutate(organism = str_replace(genome, "(\\S*)(_)(\\S*)(_)(\\S*_)(\\S*)(_)(\\S*_)(\\S*)", "\\1\\2\\3\\4\\6")) %>%
      dplyr::relocate(organism)
    epigloea_annotations <- rbindlist(list(epigloea_annotations, temp_data), use.names = T)
  }
  
epigloea_annotations_gained <- data.frame()
for (i in 1:length(gained_by_epigloea)){
  temp_data <- epigloea_annotations %>%
  dplyr::filter(grepl(gained_by_epigloea[i], Orthogroup))
  epigloea_annotations_gained <- rbindlist(list(epigloea_annotations_gained, temp_data), use.names = T)
}

epigloea_annotations_faa <- epigloea_annotations_gained %>%
  select(Transcript_ID, Translation)
colnames(epigloea_annotations_faa) = c("seq.name", "seq.text")

epigloea_gained <- dat2fasta(epigloea_annotations_faa,outfile = "gained_and_lost_strict/epigloea_orthogroups_gained.faa")
write.table(epigloea_annotations_gained, "gained_and_lost_strict/epigloea_orthogroups_gained_annotations.txt", row.names = FALSE, sep = "\t", col.names = TRUE, quote = FALSE)

others_list = c("Sporothrix_bragantina_S_CBS47491_GCA_XXXXXXXXX",
                "Sporothrix_brasiliensis_S_5110_GCA_000820605",
                "Sporothrix_curviconia_S_CBS95973_GCA_XXXXXXXXX",
                "Sporothrix_dimorphospora_S_CBS55374_GCA_021397985",
                "Sporothrix_eucalyptigena_S_CBS139899_GCA_XXXXXXXXX",
                "Sporothrix_eucalyptigena_S_CBS140593_GCA_XXXXXXXXX",
                "Sporothrix_euskadiensis_S_VPRI43754_GCA_019925375",
                "Sporothrix_globosa_S_CBS120340_GCA_001630435",
                "Sporothrix_humicola_S_CBS118129_GCA_021396245",
                "Sporothrix_inflata_S_CBS23968_GCA_021396225",
                "Sporothrix_luriei_S_CBS93772_GCA_021398005",
                "Sporothrix_mexicana_S_CBS120341_GCA_021396375",
                "Sporothrix_nigrograna_S_VPRI43755_GCA_019925305",
                "Sporothrix_pallida_S_CBS13156_GCA_021396235",
                "Sporothrix_phasma_S_CBS119721_GCA_011037845",
                "Sporothrix_protearum_S_CBS116654_GCA_016097115",
                "Sporothrix_pseudoabietina_S_VPRI43531_GCA_019925295",
                "Sporothrix_schenckii_S_1099_GCA_000961545",
                "Sporothrix_thermara_S_CBS139747MAG_GCA_XXXXXXXXX",
                "Sporothrix_variecibatus_S_CBS121961_GCA_021396255")

others_annotations <- data.frame()
for (i in 1:length(others_list)){
  temp_data <- read.delim(paste("~/Documents/2023_02_Sporothrix/10_annotations/annotation_tables_merged/",others_list[i],".annotations.merged.txt", sep = ""), stringsAsFactors = F)
  temp_data <- temp_data %>%
    #select(Transcript_ID,Orthogroup,Translation) %>%
    add_column(genome = others_list[i], .before = "Transcript_ID") %>%
    dplyr::mutate(genome = str_replace(genome, "(\\S*)(.annotations.txt)", "\\1")) %>%
    dplyr::mutate(organism = str_replace(genome, "(\\S*)(_)(\\S*)(_)(\\S*_)(\\S*)(_)(\\S*_)(\\S*)", "\\1\\2\\3\\4\\6")) %>%
    dplyr::relocate(organism)
  others_annotations <- rbindlist(list(others_annotations, temp_data), use.names = T)
}

epigloea_annotations_lost <- data.frame()
for (i in 1:length(lost_by_epigloea)){
  temp_data <- others_annotations %>%
    dplyr::filter(grepl(lost_by_epigloea[i], Orthogroup))
  epigloea_annotations_lost <- rbindlist(list(epigloea_annotations_lost, temp_data), use.names = T)
}

others_annotations_faa <- epigloea_annotations_lost %>%
  select(Transcript_ID, Translation)
colnames(others_annotations_faa) = c("seq.name", "seq.text")
epigloea_lost <- dat2fasta(others_annotations_faa,outfile = "gained_and_lost_strict/epigloea_orthogroups_lost.faa")
write.table(epigloea_annotations_lost, "gained_and_lost_strict/epigloea_orthogroups_lost_annotations.txt", row.names = FALSE, sep = "\t", col.names = TRUE, quote = FALSE)

lost_by_epigloea


################## KinFin results ##########


#load results

GO_rep <- read.delim("../KinFin/GO.cluster_functional_annotation.all.p30.x2.tsv", sep = "\t", stringsAsFactors = F, col.names = c("cluster_id", "protein_count", "taxon_count", "GO_domain_ids", "GO_domain_description"))
IPR_rep <- read.delim("../KinFin/IPR.cluster_functional_annotation.all.p30.x2.tsv", sep = "\t", stringsAsFactors = F, col.names = c("cluster_id", "protein_count", "taxon_count", "IPR_domain_ids", "IPR_domain_description"))
Pfam_rep <- read.delim("../KinFin/Pfam.cluster_functional_annotation.all.p30.x2.tsv", sep = "\t", stringsAsFactors = F, col.names = c("cluster_id", "protein_count", "taxon_count", "Pfam_domain_ids", "Pfam_domain_description"))

#combine the 3 annotation methods
orthogroup_ann <- left_join(GO_rep, IPR_rep)
orthogroup_ann <- left_join(orthogroup_ann, Pfam_rep)

write.table(orthogroup_ann, "../KinFin/combined_cluster_functional_annotation.all.p30.x2.tsv", row.names = FALSE, sep = "\t", col.names = TRUE, quote = FALSE)


#orthogroup annotations of those lost and gained by epigloea

orthogroup_ann_lost <- orthogroup_ann[orthogroup_ann$cluster_id %in% lost_by_epigloea,]

#write.table(orthogroup_ann_lost, "../KinFin/combined_cluster_functional_annotation.lost.p30.x2.tsv", row.names = FALSE, sep = "\t", col.names = TRUE, quote = FALSE)


orthogroup_ann_gained <- orthogroup_ann[orthogroup_ann$cluster_id %in% gained_by_epigloea,]
#write.table(orthogroup_ann_gained, "../KinFin/combined_cluster_functional_annotation.gained.p30.x2.tsv", row.names = FALSE, sep = "\t", col.names = TRUE, quote = FALSE)

########### orthogroup annotations were categorized manually.  Interpro ########

library(readxl)
orthogroup_ann_lost <- read_excel("../KinFin/combined_cluster_functional_annotation.lost.p30.x2.xls")

orthogroup_ann_lost_count <- orthogroup_ann_lost %>%
  count(category2) %>%
  arrange(n) %>%
  mutate(category2=factor(category2, levels=category2))

orthogroup_ann_gained <- read_excel("../KinFin/combined_cluster_functional_annotation.gained.p30.x2.xls")
orthogroup_ann_gained_count <- orthogroup_ann_gained %>%
  count(category2) %>%
  arrange(n) %>%
  mutate(category2=factor(category2, levels=category2))

orthogroup_ann_count <- orthogroup_ann_lost_count %>% 
  full_join(orthogroup_ann_gained_count, by="category2", suffix = c("_lost", "_gained")) %>%
  arrange(n_lost) %>%
  mutate(category2=factor(category2, levels=category2)) %>%
  replace_na(list(n_gained = 0))

plot_lollipop <- function(data, x, y) {
  ggplot(data, aes(x={{x}}, y={{y}})) +
    geom_segment(aes(x={{x}}, xend={{x}}, y=0, yend={{y}}), color="grey") +
    geom_point( color="orange", size=3) +
    theme_light() +
    coord_flip() +
    #theme_bw() +
    xlab("") +
    ylab("number of orthogroups") +
    theme(
      panel.grid.major.y = element_blank(),
      panel.border = element_blank(),
      axis.ticks.y = element_blank()
    )
}

p1 <- plot_lollipop(orthogroup_ann_gained_count, category2, n)
p2 <- plot_lollipop(orthogroup_ann_lost_count, category2, n)

library(gridExtra)
grid.arrange(p1, p2, nrow = 1)


p3 <- ggplot(orthogroup_ann_count) +
  geom_segment(aes(x=category2, xend=category2, y=n_lost, yend=n_gained), color="grey") +
  geom_point(aes(x=category2, y=n_lost), color="orange", size=3) +
  geom_point(aes(x=category2, y=n_gained), color="blue", size=3) +
  coord_flip() +
  theme_light() +
  theme(
    legend.position = "right",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank()) +
  xlab("") +
  ylab("number of orthogroups")


pdf(file="~/Documents/2023_02_Sporothrix/results/figures/orthogroup_annotations.pdf",
    width=5, height=4)

p3

dev.off()

str(orthogroup_ann_count_tidy)


########### orthogroup annotations were categorized manually.  GO terms ########

library(readxl)
orthogroup_ann_lost <- read_excel("../KinFin/combined_cluster_functional_annotation.lost.p30.x2.xls")

orthogroup_ann_lost_count <- orthogroup_ann_lost %>%
  separate_longer_delim(GO_domain_description, delim = ";") %>%
  count(GO_domain_description) %>%
  arrange(n) %>%
  mutate(GO_domain_description=factor(GO_domain_description, levels=GO_domain_description))

orthogroup_ann_lost_count_short <- orthogroup_ann_lost_count %>%
  filter(n > 1)

orthogroup_ann_gained <- read_excel("../KinFin/combined_cluster_functional_annotation.gained.p30.x2.xls")
orthogroup_ann_gained_count <- orthogroup_ann_gained %>%
  separate_longer_delim(GO_domain_description, delim = ";") %>%
  count(GO_domain_description) %>%
  arrange(n) %>%
  mutate(GO_domain_description=factor(GO_domain_description, levels=GO_domain_description))

orthogroup_ann_count <- orthogroup_ann_lost_count %>% 
  full_join(orthogroup_ann_gained_count, by="GO_domain_description", suffix = c("_lost", "_gained")) %>%
  replace_na(list(n_gained = 0, n_lost = 0)) %>%
  arrange(n_lost) %>%
  mutate(GO_domain_description=factor(GO_domain_description, levels=GO_domain_description))

plot_lollipop <- function(data, x, y) {
  ggplot(data, aes(x={{x}}, y={{y}})) +
    geom_segment(aes(x={{x}}, xend={{x}}, y=0, yend={{y}}), color="grey") +
    geom_point( color="orange", size=3) +
    theme_light() +
    coord_flip() +
    #theme_bw() +
    xlab("") +
    ylab("number of orthogroups") +
    theme(
      panel.grid.major.y = element_blank(),
      panel.border = element_blank(),
      axis.ticks.y = element_blank()
    )
}

p1 <- plot_lollipop(orthogroup_ann_gained_count, GO_domain_description, n)
p2 <- plot_lollipop(orthogroup_ann_lost_count, GO_domain_description, n)
p3 <- plot_lollipop(orthogroup_ann_lost_count_short, GO_domain_description, n)

library(gridExtra)

pdf(file="~/Documents/2023_02_Sporothrix/results/figures/orthogroup_annotations_GO_separate.pdf",
    width=15, height=10)

grid.arrange(p1, p3, nrow = 1)

dev.off()


p4 <- ggplot(orthogroup_ann_count) +
  geom_segment(aes(x=GO_domain_description, xend=GO_domain_description, y=n_lost, yend=n_gained), color="grey") +
  geom_point(aes(x=GO_domain_description, y=n_lost), color="orange", size=3) +
  geom_point(aes(x=GO_domain_description, y=n_gained), color="blue", size=3) +
  coord_flip() +
  theme_light() +
  theme(
    legend.position = "right",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank()) +
  xlab("") +
  ylab("number of orthogroups")


pdf(file="~/Documents/2023_02_Sporothrix/results/figures/orthogroup_annotations_GO.pdf",
    width=10, height=10)

p4

dev.off()


