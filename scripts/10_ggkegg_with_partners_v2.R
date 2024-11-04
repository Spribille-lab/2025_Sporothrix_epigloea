#BiocManager::install("ggkegg")
#devtools::install_github("noriakis/ggkegg")
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
setwd("~/Documents/2023_02_Sporothrix/10_annotations/kegg_relaxed/")

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
write.table(hdf_combined, "../../results/tables/modulecompleteness_relaxed_combined.txt", sep="\t", quote=FALSE, row.names = TRUE, col.names = NA)

str(hdf_combined)

warnings()
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
pdf(file="~/Documents/2023_02_Sporothrix/results/figures/module_completeness_relaxed_with_partners.pdf",
    width=13, height=30)
draw(ht2, heatmap_legend_side = "left")
dev.off()






## 1. misc
library(tidyverse)
library(stringr)

annotated_genomes<-read.delim("../genomes_with_partners.txt",header=T,col.names = c("Genome", "Group"))

##### 1. Combine kegg mapper files for all genomes#####
read_kegg_table<-function(genome_list){
  
  read_kegg_by_genome<-function(genome){
    filename<-paste0("kegg_output/",genome,"/",genome,".kegg.mapper.txt")
    keggfile<-read.delim(filename,header=T,na.strings = "", col.names = c("GeneID", "KO"))
    keggfile$Genome<-genome
    return(keggfile)}
  
  l<-lapply(genome_list,read_kegg_by_genome)
  
  return(l)   
}

l<-read_kegg_table(annotated_genomes$Genome)
warnings()
kegg_combined<-do.call(rbind,l)
#write.table(kegg_combined, "summarized_outputs/kegg_combined.txt", sep='\t',quote = F, row.names = F, col.names = T)


#2. kegg annotations
###load data
kegg_of_interest<-read.delim("kegg_of_interest_fig.txt")
#kegg_combined<-read.delim("summarized_outputs/kegg_combined.txt")
colnames(kegg_combined)<-c("locustag","ID","Genome")



###prepare dataset
kegg_list<-kegg_of_interest[,1]

kegg_df <- kegg_combined %>% filter(ID %in% kegg_list) %>% group_by(Genome,ID) %>%
  summarize(sum=n()) %>% mutate(completeness=ifelse(sum>0,1,0)) %>% dplyr::select(-sum) %>%
  pivot_wider(names_from = ID,values_from = completeness, values_fill = 0) %>%
  select(Genome, K02429, K03320) %>%
  #new variable: Fucose transporter
  #mutate(fucose_transporter = ifelse(K02429==1,1,0)) %>%
  #new variable: ammonium transporter
  #mutate(ammonium_transporter = ifelse(K03320==1,1,0)) %>%
  #remove ID columns
  #dplyr::select(-contains("K"))# %>% left_join( annotated_mags_arranged)
  pivot_longer(-Genome, names_to = "ID", values_to = "completeness") %>%
  mutate(pathway = case_when(ID == "K02429" ~ "fucose_transporter",
                             ID == "K03320" ~ "ammonium_transporter"))
#############
# 5. vitamin modules
module_kegg<-read.delim("../../results/tables/modulecompleteness_relaxed_combined.txt", header=T)

module_kegg_long <- module_kegg %>% pivot_longer(-X, names_to = "Genome",values_to = "completeness") %>% 
  mutate(module=gsub("\\_.*","",X))
colnames(module_kegg_long)[1] <- "ID"

modules_kegg_sel <- module_kegg_long %>% filter(module %in% c("M00899","M00898",
                                                              "M00572", "M00123",
                                                              "M00911")) %>%
  group_by(Genome,module) %>%
  #mutate(mean_completeness=mean(completeness)) %>%
  mutate(pathway=case_when(module == "M00898" ~"thiamine_synthesis",
                           module == "M00899" ~ "thiamine_salvage",
                           module %in% c("M00572","M00123") ~ "biotin",
                           module == "M00911" ~ "riboflavin")) %>%
  #select(-X) %>%
  #distinct() %>%
  group_by(Genome,pathway) %>%
  #mutate(presence=completeness) %>%
  #mutate(presence=case_when(completeness==1 ~ "full",
  #                         completeness<1&completeness>0.9 ~ "partial",
  #                        completeness<=0.9 ~ "missing")) %>%
  #select(-c(completeness,module)) %>%
  select(-module) %>%
  distinct()

# 6.  viz

#df_long <- kegg_df %>%
#  rbind(modules_kegg_sel)

df_long <- modules_kegg_sel

df_long <- annotated_genomes %>% 
  inner_join(df_long) %>% mutate(presence_factor = ifelse(
    completeness==1,"full",ifelse(completeness>0.5,"partial","missing"))) %>%
  mutate(presence_size = ifelse(
    completeness>0,1,NA)) #%>%
#mutate(bac_family_label = ifelse(bac_family=="Sphingomonadaceae","Sphingo\nmonadaceae",
#                                 ifelse(bac_family=="UBA10450","UBA\n10450",
#                                        ifelse(bac_family=="Beijerinckiaceae","Beijerin\nckiaceae",bac_family))))

#order rows and columns
df_long$Genome<-factor(df_long$Genome,level = annotated_genomes$Genome)
str(df_long)

df_long$pathway %>% unique

pathway_type<-data.frame("pathway"=c("capsular_transporter",
                                     "branched_transporter" ,         
                                     "l_amino_transporter",           
                                     "glutamate_transporter" ,
                                     "ammonium_transporter",
                                     "urea_transporter",
                                     "ribose_transporter" ,           
                                     "xylose_transporter" ,           
                                     "multiple_sugar_transporter"  ,  
                                     "fructose_transporter" ,         
                                     "arabinose_transporter" ,
                                     "fucose_transporter",
                                     "erythritol_transporter",        
                                     "xylitol_transporter",           
                                     "inositol_transporter",          
                                     "glycerol_transporter",
                                     "glycerol_aquaporin_transporter",
                                     "glycerol_sorbitol_transporter",                      
                                     "sorbitol_mannitol_transporter" ,              
                                     "thiamine_synthesis",
                                     "thiamine_salvage",
                                     "riboflavin"  ,   
                                     "cobalamin",
                                     "biotin" ,
                                     "urease",
                                     "nitrogen_fixation",
                                     "methanol_dehydrogenase",
                                     "Calvin_cycle",                  
                                     "carotenoids",
                                     "bacteriochlorophyll",
                                     "photosystem"
),"functions"=c(rep("Other transporters",6),rep("Carbohydrate transporters",13),
                rep("Thiamine\nsynthesis",1), rep("Thiamine\nsalvage",1), rep("Riboflavin\nsynthesis",1), 	
                rep("Cobalamin", 1),
                rep("Biotin\nsynthesis", 1),
                rep("C and N metabolism",4), rep("Photo-\nsynthesis",3)
))

df_long <- df_long %>% left_join(pathway_type)
df_long$pathway <- factor(df_long$pathway, level = c("biotin",
                                                     "thiamine_synthesis",
                                                     "thiamine_salvage",
                                                     "riboflavin",
                                                     "capsular_transporter",
                                                     "branched_transporter" ,         
                                                     "l_amino_transporter",           
                                                     "glutamate_transporter" ,
                                                     "ammonium_transporter",
                                                     "urea_transporter",
                                                     "ribose_transporter" ,           
                                                     "xylose_transporter" ,           
                                                     "multiple_sugar_transporter"  ,  
                                                     "fructose_transporter" ,         
                                                     "arabinose_transporter" ,
                                                     "fucose_transporter",
                                                     "erythritol_transporter",        
                                                     "xylitol_transporter",           
                                                     "inositol_transporter",          
                                                     "glycerol_transporter",
                                                     "glycerol_aquaporin_transporter",
                                                     "glycerol_sorbitol_transporter",                      
                                                     "sorbitol_mannitol_transporter",
                                                     "cobalamin",
                                                     "iron_ion_transport","siderophore_synthesis",
                                                     "urease",
                                                     "nitrogen_fixation",
                                                     "methanol_dehydrogenase",
                                                     "Calvin_cycle",                  
                                                     "carotenoids",
                                                     "bacteriochlorophyll",
                                                     "photosystem",
                                                     "exopolysaccaride"
))

df_long$functions <- factor(df_long$functions,level = c("Photo-\nsynthesis","C and N metabolism","Biotin\nsynthesis","Thiamine\nsynthesis","Thiamine\nsalvage","Riboflavin\nsynthesis","Cobalomin","Carbohydrate transporters", "Other transporters"))
df_long$Group <- factor(df_long$Group,level = c("core Sporothrix clade", "outgroup", "associates"))
df_long <- df_long %>% arrange(ID)
df_long$ID <- factor(df_long$ID, level = rev(unique(df_long$ID) )     )
df_long$Genome %>% unique()                     
str(df_long)

cols = c("#2a778dff", "#22a985ff", "#440055ff")

pdf(file="~/Documents/2023_02_Sporothrix/results/figures/kegg_relaxed.pdf",
    width=18, height=8)

ggplot(df_long, aes(x=Genome,y=ID,size=presence_size,shape=presence_factor,color=Group)) +
  geom_point(size=4) +
  scale_colour_manual(values = cols) +
  scale_shape_manual(values=c(16,10), limits = c("full","partial")) +
  guides(size="none",color="none",shape = "none") + theme_minimal()+
  facet_grid(functions~Group,scales = "free", space="free", labeller = label_wrap_gen(width=5)) +
  theme(#axis.text.x = element_blank(),
    axis.text.x = element_text(size =12, angle = 90, hjust=1),
    axis.text.y = element_text(size=12),
    strip.text.x = element_text(size =12),
    strip.text.y = element_text(size =12, angle = 0),
    panel.spacing = unit(2, "lines"),
    legend.title = element_text(size=7),
    legend.text = element_text(size=7)) +
  xlab("") +
  ylab("") +
  scale_x_discrete(labels=c("Sporothrix_euskadiensis_S_VPRI43754_GCA_019925375" = "S. euskadiensis VPRI43754",
                            "Sporothrix_pseudoabietina_S_VPRI43531_GCA_019925295" = "S. pseudoabietina VPRI43531",
                            "Sporothrix_variecibatus_S_CBS121961_GCA_021396255" = "S. variecibatus CBS121961",
                            "Sporothrix_protearum_S_CBS116654_GCA_016097115" = "S. protearum CBS116654",
                            "Sporothrix_humicola_S_CBS118129_GCA_021396245" = "S. humicola CBS118129",
                            "Sporothrix_pallida_S_CBS13156_GCA_021396235" = "S. pallida CBS13156",    
                            "Sporothrix_mexicana_S_CBS120341_GCA_021396375" = "S. mexicana CBS120341",
                            "Sporothrix_brasiliensis_S_5110_GCA_000820605" = "S. brasiliensis 5110",  
                            "Sporothrix_schenckii_S_1099_GCA_000961545" = "S. schenckii 1099",
                            "Sporothrix_globosa_S_CBS120340_GCA_001630435" = "S. globosa CBS120340",  
                            "Sporothrix_luriei_S_CBS93772_GCA_021398005" = "S. luriei CBS937.72",
                            "Sporothrix_phasma_S_CBS119721_GCA_011037845" = "S. phasma CBS119721",
                            "Sporothrix_dimorphospora_S_CBS55374_GCA_021397985" = "S. dimorphospora CBS553.74",
                            "Sporothrix_inflata_S_CBS23968_GCA_021396225" = "S. inflata CBS239.68",
                            "Sporothrix_bragantina_S_CBS47491_GCA_XXXXXXXXX" = "S. bragantina CBS47491",
                            "Sporothrix_curviconia_S_CBS95973_GCA_XXXXXXXXX" = "S. curviconia CBS95973",     
                            "Sporothrix_thermara_S_CBS139747_GCA_XXXXXXXXX" =  "S. thermara CBS139747",
                            "Sporothrix_epigloea_S_CBS57363_GCA_943900835" = "S. epigloea CBS573.63",       
                            "Sporothrix_epigloea_S_TF4163_GCA_944036445" = "S. epigloea TF4163",
                            "Sporothrix_epigloea_S_CBS119000_GCA_943908295" = "S. epigloea CBS119000",      
                            "Sporothrix_eucalyptigena_S_CBS139899_GCA_XXXXXXXXX" = "S. eucalyptigena CBS139899",
                            "Sporothrix_eucalyptigena_S_CBS140593_GCA_XXXXXXXXX" = "S. eucalyptigena CBS140593",  
                            "Sporothrix_nigrograna_S_VPRI43755_GCA_019925305" = "S. nigrograna VPRI43755",
                            "Ophiostoma_novoulmi_S_H327_GCA_000317715" = "O. novoulmi H327",      
                            "Ophiostoma_ips_S_VPRI43529_GCA_019925475" = "O. ips VPRI43529",
                            "Sporothrix_brunneoviolacea_S_CBS124561_GCA_021396205" = "S. brunneoviolacea CBS124561",
                            "Ophiostoma_fasciatum_S_VPRI43845_GCA_019925495" = "O. fasciatum VPRI43845",
                            "Sporothrix_insectorum_S_RCEF264_GCA_001636815" = "S. insectorum RCEF264",      
                            "Leptographium_lundbergii_S_CBS138716_GCA_001455505" = "L. lundbergii CBS138716",
                            "Tremella_yokohamensis_S_CBS18117_GCA_963924285" = "T. fuciformis CBS18117",
                            "Annulohypoxylon_annulatum_S_CBS149473_GCA_963924275" = "A. annulatum CBS149473")) +
  scale_y_discrete(labels=c("M00899_2" = "HMP => TMP | M00899_2",
                            "M00899_1" = "HET => TMP | M00899_1",
                            "M00898_2" = "pyridoxal phosphate => TMP/thiamine/TPP | M00898_2",
                            "M00898_1" = "NAD+ + Glycine + [Protein]-L-cysteine => HET | M00898_1",
                            "M00911_3" = "intermediates => riboflavin/FMN/FAD | M00911_3",
                            "M00911_2" = "ribulose 5-P => dihydroxybutan-2-one 4-phosphate | M00911_2",
                            "M00911_1" = "GTP => 5-amino-6-(ribitylamino)uracil | M00911_1",
                            "M00572_1" = "BioC-BioH pathway, malonyl-ACP => pimeloyl-ACP | M00572",
                            "M00123_1" = "pimeloyl-ACP => biotin | M00123"))

dev.off()

# module annotations
moddesc <- data.table::fread("https://rest.kegg.jp/list/module", header=FALSE)
colnames(moddesc) = c("Module", "Description")
modrank <- data.table::fread("module_ranks.txt", header=TRUE) # how the modules are organized at https://www.genome.jp/brite/ko00002
moddesrank <- dplyr::left_join(modrank, moddesc, by = "Module")


# which genomes had BioF
BioF <- kegg_combined %>% select(ID, Genome) %>% filter(ID != "NA") %>% filter(ID == "K00652") %>% unique()
BioF <- left_join(annotated_genomes, BioF)
# missing in protearum

# which genomes had Thi5
Thi5 <- kegg_combined %>% select(ID, Genome) %>% filter(ID != "NA") %>% filter(ID == "K18278") %>% unique()
Thi5 <- left_join(annotated_genomes, Thi5)

# which genomes had ThiD
ThiD <- kegg_combined %>% select(ID, Genome) %>% filter(ID != "NA") %>% filter(ID == "K00941") %>% unique()
ThiD <- left_join(annotated_genomes, ThiD)

# which genomes had Rib2
Rib2 <- kegg_combined %>% select(ID, Genome) %>% filter(ID != "NA") %>% filter(ID == "K14655") %>% unique()
Rib2 <- left_join(annotated_genomes, Rib2)

