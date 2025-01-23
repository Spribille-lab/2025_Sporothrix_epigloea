
# make an ultrametric tree for CAFE analysis

setwd("~/Documents/2023_02_Sporothrix/14_CAFE/")
library(ape)
library(phytools)
library(tidyverse)
library(ggtree)
library(treeio)
library(tidytree)

ML_tree <- ape::read.tree("~/Documents/2023_02_Sporothrix/09_orthology/OrthoFinder/Results_Aug03/Species_Tree/iqtree2/sporothrix_placement.rerooted.contree")
class(ML_tree)
ML_tree$tip.label
ML_tree <- drop.tip(ML_tree, 25)
str(ML_tree)


ML_tree <- ape::chronos(ML_tree, lambda=1)


#Setting initial dates...
#Fitting in progress... get a first set of estimates
#(Penalised) log-lik = -10.12795 
#Optimising rates... dates... -10.12795 
#Optimising rates... dates... -10.12784 

#log-Lik = -10.10604 
#PHIIC = 186.44 

is.binary(ML_tree)
is.ultrametric(ML_tree)
is.rooted(ML_tree)

plotTree(ML_tree,fsize=1) 
#write.tree(ML_tree, file="~/Documents/2023_02_Sporothrix/09_orthology/OrthoFinder/Results_Aug03/Species_Tree/iqtree2/sporothrix_placement.rerooted.ultrametric.contree")

# format protease family results from funannotate for CAFE5

protease_family <- read.delim("~/Documents/2023_02_Sporothrix/10_annotations/funannotate/funannotate_compare_20230809/merops/MEROPS.all.results.csv", header = T, na.strings='', sep=",") %>%
  mutate(Desc = "(null)") %>%
  relocate(Desc)
colnames(protease_family)[2] <- "Family ID"
#write.table(protease_family, file="~/Documents/2023_02_Sporothrix/14_CAFE/protease_families_funannotate.txt", sep='\t', quote = F, row.names = F, col.names = T)

# format protease class results from funannotate for CAFE5
# in the end I did not use this

protease_class <- read.delim("~/Documents/2023_02_Sporothrix/10_annotations/funannotate/funannotate_compare_20230809/merops/MEROPS.summary.results.csv", header = T, na.strings='', sep=",") %>%
  mutate(Desc = "(null)") %>%
  relocate(Desc)
colnames(protease_class)[2] <- "Family ID"
#write.table(protease_class, file="~/Documents/2023_02_Sporothrix/14_CAFE/protease_classes_funannotate.txt", sep='\t', quote = F, row.names = F, col.names = T)

# format secmets class results from antiSMASH for CAFE5
secmets<-read.delim("~/Documents/2023_02_Sporothrix/10_annotations/antiSMASH/antismash.txt") %>% 
  dplyr::select(genome, secondary_metabolite_type) %>%
  dplyr::group_by(genome) %>%
  dplyr::count(secondary_metabolite_type) %>%
  tidyr::spread(secondary_metabolite_type, n, fill=0) %>%
  base::as.data.frame()
row.names(secmets) = secmets$genome
secmets <- secmets %>% dplyr::select(-genome)
secmets <- t(secmets)
secmets <- as.data.frame(secmets)
secmets <- secmets %>% rownames_to_column(var = "Family ID") %>%
  mutate(Desc = "(null)") %>%
  relocate(Desc)
#write.table(secmets, file="~/Documents/2023_02_Sporothrix/14_CAFE/secmets.txt", sep='\t', quote = F, row.names = F, col.names = T)

# make tree for orthogroups and format orthofinder output for CAFE

# remove outgroups that were not included in the orthogroup inference
ML_tree$tip.label
ML_tree_orthog <- drop.tip(ML_tree, c(24:29))
plotTree(ML_tree_orthog,fsize=1)
ggtree(ML_tree) + geom_text(aes(label=node), hjust=-.3) + geom_tiplab(hjust=0.2)
ggtree(ML_tree_orthog) + geom_text(aes(label=node), hjust=-.3) + geom_tiplab(hjust=1)
is.binary(ML_tree_orthog)
is.ultrametric(ML_tree_orthog)
is.rooted(ML_tree_orthog)
#write.tree(ML_tree_orthog, file="~/Documents/2023_02_Sporothrix/09_orthology/OrthoFinder/Results_Aug03/Species_Tree/iqtree2/sporothrix_placement.rerooted.ultrametric.coresporothrix.contree")
ML_tree_orthog <- read.tree("~/Documents/2023_02_Sporothrix/09_orthology/OrthoFinder/Results_Aug03/Species_Tree/iqtree2/sporothrix_placement.rerooted.ultrametric.coresporothrix.contree")
ggtree(ML_tree_orthog) + geom_text(aes(label=node), hjust=-.3) + geom_tiplab(hjust=1)

orthogroups <- read.delim("~/Documents/2023_02_Sporothrix/09_orthology/OrthoFinder/Results_Aug07/Orthogroups/Orthogroups.GeneCount.tsv", header = T, na.strings='', sep="\t") %>%
  mutate(Desc = "(null)") %>%
  relocate(Desc) %>%
  select(-Total)
colnames(orthogroups)[2] <- "Family ID"
#write.table(orthogroups, file="~/Documents/2023_02_Sporothrix/14_CAFE/orthogroups.txt", sep='\t', quote = F, row.names = F, col.names = T)

# visualize CAFE results





tree <- read.tree("~/Documents/2023_02_Sporothrix/09_orthology/OrthoFinder/Results_Aug03/Species_Tree/iqtree2/sporothrix_placement.rerooted.ultrametric.contree")

#draw tree with node numbers
ggtree(tree) + geom_text(aes(label=node), hjust=-.3) + geom_tiplab(hjust=-0.3)

# convert tree into tibble tree
tbl_tree <- as_tibble(tree)
as.treedata(tbl_tree)
#ggtree(tbl_tree)
str(tbl_tree)

# extract number of families from CAFE analysis
node_data <- read.delim("orthogroups_singlelambda_k1/Base_clade_results.txt", header=T) %>%
  as.data.frame() %>%
  mutate(node = str_replace(X.Taxon_ID, "(^.*<)(\\d*)(>)","\\2")) %>%
  select(-X.Taxon_ID) %>%
  relocate(node)
node_data$node <- as.integer(node_data$node)
node_data <- node_data  %>%
  arrange(node)
str(node_data)

# join with tree
tree_data <- full_join(tbl_tree, node_data, by = 'node')

# make into tree data object
tree_data <- as.treedata(tree_data)
str(tree_data)

tree_labels <- data.frame(label = ML_tree$tip.label,
                          label2 = c("S. euskadiensis VPRI43754",
                                     "S. pseudoabietina VPRI43531",
                                     "S. variecibatus CBS 121961",
                                     "S. protearum CBS 116654",
                                     "S. humicola CBS 118129",
                                     "S. pallida CBS 131.56",
                                     "S. mexicana CBS 120341",
                                     "S. brasiliensis 5110",
                                     "S. schenckii 1099-18",
                                     "S. globosa CBS 120340",
                                     "S. luriei CBS 937.72",
                                     "S. phasma CBS 119721",
                                     "S. dimorphospora CBS 553.74",
                                     "S. inflata CBS 239.68",
                                     "S. bragantina CBS 474.91",
                                     "S. curviconia CBS 959.73",
                                     "S. thermara CBS 139747 MAG",
                                     "S. epigloea CBS 573.63",
                                     "S. epigloea TF4-1_63 MAG",
                                     "S. epigloea CBS 119000",
                                     "S. eucalyptigena CBS 139899",
                                     "S. eucalyptigena CBS 140593",
                                     "S. nigrograna VPRI43755",
                                     "O. novo-ulmi H327",
                                     "O. ips VPRI43529",
                                     "S. brunneoviolacea CBS 124561",
                                     "O. fasciatum VPRI43845",
                                     "S. insectorum RCEF264",
                                     "L. lundbergii CBS 138716"))

tree_data <- left_join(tree_data, tree_labels, by = 'label')
#child(tree_data, 51)

line_type_full_tree <- data.frame(node=1:Nnode2(ML_tree), lty = 1)
line_type_full_tree[c(54,55,18,19,20), 2] <- c(4,4,4,4,4)

line_type_core_tree <- data.frame(node=1:Nnode2(ML_tree_orthog), lty = 1)
line_type_core_tree[c(43,44,18,19,20), 2] <- c(4,4,4,4,4)

ggtree(tree_data) %<+% line_type_full_tree + aes(linetype=I(lty)) +
  geom_tiplab(aes(label=label2),offset=.05) +
  theme_tree2() +
  ggplot2::xlim(0, 2) +
  geom_point(aes(size = Decrease)) +
  scale_size(range = c(0, 11))





# use the above to create plotting functions

plot_CAFE_decrease <- function(target_tree, target_data, target_group, xlim, target_line_type){

tree <- read.tree(target_tree)

#draw tree with node numbers
#ggtree(tree) + geom_text(aes(label=node), hjust=-.3) + geom_tiplab(hjust=-0.3)

# convert tree into tibble tree
tbl_tree <- as_tibble(tree)
as.treedata(tbl_tree)

# extract number of families from CAFE analysis
node_data <- read.delim(target_data, header=T) %>%
  as.data.frame() %>%
  mutate(node = str_replace(X.Taxon_ID, "(^.*<)(\\d*)(>)","\\2")) %>%
  select(-X.Taxon_ID) %>%
  relocate(node)
node_data$node <- as.integer(node_data$node)
node_data <- node_data  %>%
  arrange(node)
str(node_data)

# join with tree
tree_data <- full_join(tbl_tree, node_data, by = 'node')

# make into tree data object
tree_data <- as.treedata(tree_data)
str(tree_data)

tree_data <- left_join(tree_data, tree_labels, by = 'label')

p <- ggtree(tree_data) %<+% target_line_type + aes(linetype=I(lty)) +
  geom_tiplab(aes(label=label2), offset=.05, size = 3) +
  theme_tree() +
  ggplot2::xlim(0, xlim) +
  vexpand(.05, direction = -1) +
  geom_point(aes(size = Decrease), colour="skyblue4", alpha=0.5) +
  scale_size(range = c(0, 11)) +
  ggtitle(paste("Contractions - ",target_group,sep = "")) +
  theme(plot.title = element_text(hjust = 0.5))

return(p)

}

plot_CAFE_increase <- function(target_tree, target_data, target_group, xlim, target_line_type){
  
  tree <- read.tree(target_tree)
  
  #draw tree with node numbers
  #ggtree(tree) + geom_text(aes(label=node), hjust=-.3) + geom_tiplab(hjust=-0.3)
  
  # convert tree into tibble tree
  tbl_tree <- as_tibble(tree)
  as.treedata(tbl_tree)
  
  # extract number of families from CAFE analysis
  node_data <- read.delim(target_data, header=T) %>%
    as.data.frame() %>%
    mutate(node = str_replace(X.Taxon_ID, "(^.*<)(\\d*)(>)","\\2")) %>%
    select(-X.Taxon_ID) %>%
    relocate(node)
  node_data$node <- as.integer(node_data$node)
  node_data <- node_data  %>%
    arrange(node)
  str(node_data)
  
  # join with tree
  tree_data <- full_join(tbl_tree, node_data, by = 'node')
  
  # make into tree data object
  tree_data <- as.treedata(tree_data)
  str(tree_data)
  
  tree_data <- left_join(tree_data, tree_labels, by = 'label')
  
  p <- ggtree(tree_data) %<+% target_line_type + aes(linetype=I(lty)) +
    geom_tiplab(aes(label=label2), offset=.05, size = 3)  +
    theme_tree() +
    ggplot2::xlim(0, xlim) +
    vexpand(.05, direction = -1) +
    geom_point(aes(size = Increase), colour="antiquewhite4", alpha=0.5) +
    scale_size(range = c(0, 11)) +
    ggtitle(paste("Expansions - ",target_group,sep = "")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
  
}

cazyme_decrease <- plot_CAFE_decrease("~/Documents/2023_02_Sporothrix/09_orthology/OrthoFinder/Results_Aug03/Species_Tree/iqtree2/sporothrix_placement.rerooted.ultrametric.contree",
                                      "cazymefamily_singlelambda_k1/Base_clade_results.txt",
                                      target_group = "CAZyme families",
                                      xlim = 2,
                                      line_type_full_tree)

protease_decrease <- plot_CAFE_decrease("~/Documents/2023_02_Sporothrix/09_orthology/OrthoFinder/Results_Aug03/Species_Tree/iqtree2/sporothrix_placement.rerooted.ultrametric.contree",
                                        "proteasefamily_singlelambda_k1/Base_clade_results.txt",
                                        target_group = "Protease families",
                                        xlim = 2,
                                        line_type_full_tree)

orthogroups_decrease <- plot_CAFE_decrease("~/Documents/2023_02_Sporothrix/09_orthology/OrthoFinder/Results_Aug03/Species_Tree/iqtree2/sporothrix_placement.rerooted.ultrametric.coresporothrix.contree",
                                       "orthogroups_singlelambda_k1/Base_clade_results.txt",
                                       target_group = "Orthogroups",
                                       xlim = 1,
                                       line_type_core_tree)

cazyme_increase <- plot_CAFE_increase("~/Documents/2023_02_Sporothrix/09_orthology/OrthoFinder/Results_Aug03/Species_Tree/iqtree2/sporothrix_placement.rerooted.ultrametric.contree",
                                      "cazymefamily_singlelambda_k1/Base_clade_results.txt",
                                      target_group = "CAZyme families",
                                      xlim = 2,
                                      line_type_full_tree)

protease_increase <- plot_CAFE_increase("~/Documents/2023_02_Sporothrix/09_orthology/OrthoFinder/Results_Aug03/Species_Tree/iqtree2/sporothrix_placement.rerooted.ultrametric.contree",
                                        "proteasefamily_singlelambda_k1/Base_clade_results.txt",
                                        target_group = "Protease families",
                                        xlim = 2,
                                        line_type_full_tree)


orthogroups_increase <- plot_CAFE_increase("~/Documents/2023_02_Sporothrix/09_orthology/OrthoFinder/Results_Aug03/Species_Tree/iqtree2/sporothrix_placement.rerooted.ultrametric.coresporothrix.contree",
                                           "orthogroups_singlelambda_k1/Base_clade_results.txt",
                                           target_group = "Orthogroups",
                                           xlim = 1,
                                           line_type_core_tree)


library(cowplot)

pdf(file="~/Documents/2023_02_Sporothrix/results/figures/families_expansion_contraction.pdf", width=10, height=10)
plot_grid(cazyme_decrease, cazyme_increase, protease_decrease, protease_increase, ncol = 2)
dev.off()


pdf(file="~/Documents/2023_02_Sporothrix/results/figures/orthogroups_expansion_contraction.pdf", width=9, height=5)
plot_grid(orthogroups_decrease, orthogroups_increase, ncol = 2)
dev.off()






# process orthogroup results

# load kinfin results
GO_rep <- read.delim("../09_orthology/OrthoFinder/KinFin/GO.cluster_functional_annotation.all.p30.x2.tsv", sep = "\t", stringsAsFactors = F, col.names = c("cluster_id", "protein_count", "taxon_count", "GO_domain_ids", "GO_domain_description"))
IPR_rep <- read.delim("../09_orthology/OrthoFinder/KinFin/IPR.cluster_functional_annotation.all.p30.x2.tsv", sep = "\t", stringsAsFactors = F, col.names = c("cluster_id", "protein_count", "taxon_count", "IPR_domain_ids", "IPR_domain_description"))
Pfam_rep <- read.delim("../09_orthology/OrthoFinder/KinFin/Pfam.cluster_functional_annotation.all.p30.x2.tsv", sep = "\t", stringsAsFactors = F, col.names = c("cluster_id", "protein_count", "taxon_count", "Pfam_domain_ids", "Pfam_domain_description"))

#combine the 3 annotation methods
orthogroup_ann <- left_join(GO_rep, IPR_rep)
orthogroup_ann <- left_join(orthogroup_ann, Pfam_rep)

#####################epigloea transition################ not used in the end ####



epigloea_orthogroup_taxa = c(43,44, "Sporothrix_epigloea_CBS57363", "Sporothrix_epigloea_TF4163", "Sporothrix_epigloea_CBS119000")
epigloea_families_taxa = c(54,55, "Sporothrix_epigloea_CBS57363", "Sporothrix_epigloea_TF4163", "Sporothrix_epigloea_CBS119000")
pathogenic_orthogroup_taxa = c(34,35,36,"Sporothrix_schenckii_1099", "Sporothrix_brasiliensis_5110", "Sporothrix_globosa_CBS120340", "Sporothrix_luriei_CBS93772")
pathogenic_families_taxa = c(45,46,47,"Sporothrix_schenckii_1099", "Sporothrix_brasiliensis_5110", "Sporothrix_globosa_CBS120340", "Sporothrix_luriei_CBS93772")




orthogroup_transition <- read.delim("orthogroups_singlelambda_k1_cafeplotter/result_summary.tsv", header = T, na.strings='', sep="\t") %>%
  filter(TaxonID %in% c(43,44,34,35,36, "Sporothrix_epigloea_CBS57363", "Sporothrix_epigloea_TF4163", "Sporothrix_epigloea_CBS119000", "Sporothrix_schenckii_1099", "Sporothrix_brasiliensis_5110", "Sporothrix_globosa_CBS120340", "Sporothrix_luriei_CBS93772")) %>%
  filter(Pvalue < 0.05)

cazyme_transition <- read.delim("cazymefamily_singlelambda_k1_cafeplotter/result_summary.tsv", header = T, na.strings='', sep="\t") %>%
  filter(TaxonID %in% c(54,55,45,46,47, "Sporothrix_epigloea_CBS57363", "Sporothrix_epigloea_TF4163", "Sporothrix_epigloea_CBS119000", "Sporothrix_schenckii_1099", "Sporothrix_brasiliensis_5110", "Sporothrix_globosa_CBS120340", "Sporothrix_luriei_CBS93772")) %>%
  filter(Pvalue < 0.05) 

protease_transition <- read.delim("proteasefamily_singlelambda_k1_cafeplotter/result_summary.tsv", header = T, na.strings='', sep="\t") %>%
  filter(TaxonID %in% c(54,55,45,46,47, "Sporothrix_epigloea_CBS57363", "Sporothrix_epigloea_TF4163", "Sporothrix_epigloea_CBS119000", "Sporothrix_schenckii_1099", "Sporothrix_brasiliensis_5110", "Sporothrix_globosa_CBS120340", "Sporothrix_luriei_CBS93772")) %>%
  filter(Pvalue < 0.05) 
  





orthogroup_epigloea_transition <- read.delim("orthogroups_singlelambda_k1_cafeplotter/result_summary.tsv", header = T, na.strings='', sep="\t") %>%
  filter(TaxonID %in% c(epigloea_orthogroup_taxa)) %>%
  arrange(Pvalue) %>%
  filter(Pvalue < 0.05) %>%
  #filter(Change < 0) %>%
  left_join(orthogroup_ann, by = join_by(FamilyID == cluster_id)) %>%
  select(-protein_count, -taxon_count)
write.table(orthogroup_epigloea_transition, file="~/Documents/2023_02_Sporothrix/results/tables/epigloea_transition_orthogroups.txt", sep='\t', quote = F, row.names = F, col.names = T)

#number of rapidly evolving orthogroups that were significant for these branches
sig_og <- orthogroup_epigloea_transition %>% select(FamilyID) %>% unique()

count(orthogroup_epigloea_transition %>% select(FamilyID) %>% unique())

#number of rapidly evolving orthogroups that were significant for these branches and decreased
count(orthogroup_epigloea_transition %>% filter(Change < 0) %>% select(FamilyID) %>% unique())

#number of rapidly evolving orthogroups that were significant for these branches and increased
count(orthogroup_epigloea_transition %>% filter(Change > 0) %>% select(FamilyID) %>% unique())

count(orthogroup_epigloea_transition %>% filter(Change == 0) %>% select(FamilyID) %>% unique())

# process cazyme results

cazyme_epigloea_transition <- read.delim("cazymefamily_singlelambda_k1_cafeplotter/result_summary.tsv", header = T, na.strings='', sep="\t") %>%
  filter(TaxonID %in% c(epigloea_families_taxa)) %>%
  arrange(Pvalue) %>%
  filter(Pvalue < 0.05) %>%
  #filter(Change < 0) %>%
  left_join(orthogroup_ann, by = join_by(FamilyID == cluster_id)) %>%
  select(-protein_count, -taxon_count)
write.table(cazyme_epigloea_transition, file="~/Documents/2023_02_Sporothrix/results/tables/epigloea_transition_cazymes.txt", sep='\t', quote = F, row.names = F, col.names = T)

# process protease results


protease_epigloea_transition <- read.delim("proteasefamily_singlelambda_k1_cafeplotter/result_summary.tsv", header = T, na.strings='', sep="\t") %>%
  filter(TaxonID %in% c(epigloea_families_taxa)) %>%
  arrange(Pvalue) %>%
  filter(Pvalue < 0.05) %>%
  #filter(Change < 0) %>%
  left_join(orthogroup_ann, by = join_by(FamilyID == cluster_id)) %>%
  select(-protein_count, -taxon_count)
write.table(protease_epigloea_transition, file="~/Documents/2023_02_Sporothrix/results/tables/epigloea_transition_protease.txt", sep='\t', quote = F, row.names = F, col.names = T)

#####################pathogenic transition################


orthogroup_pathogenic_transition <- read.delim("orthogroups_singlelambda_k1_cafeplotter/result_summary.tsv", header = T, na.strings='', sep="\t") %>%
  filter(TaxonID %in% c(pathogenic_orthogroup_taxa)) %>%
  arrange(Pvalue) %>%
  filter(Pvalue < 0.05) %>%
  #filter(Change < 0) %>%
  left_join(orthogroup_ann, by = join_by(FamilyID == cluster_id)) %>%
  select(-protein_count, -taxon_count)
write.table(orthogroup_pathogenic_transition, file="~/Documents/2023_02_Sporothrix/results/tables/pathogenic_transition_orthogroups.txt", sep='\t', quote = F, row.names = F, col.names = T)

# process cazyme results

cazyme_pathogenic_transition <- read.delim("cazymefamily_singlelambda_k1_cafeplotter/result_summary.tsv", header = T, na.strings='', sep="\t") %>%
  filter(TaxonID %in% c(pathogenic_families_taxa)) %>%
  arrange(Pvalue) %>%
  filter(Pvalue < 0.05) %>%
  #filter(Change < 0) %>%
  left_join(orthogroup_ann, by = join_by(FamilyID == cluster_id)) %>%
  select(-protein_count, -taxon_count)
write.table(cazyme_pathogenic_transition, file="~/Documents/2023_02_Sporothrix/results/tables/pathogenic_transition_cazymes.txt", sep='\t', quote = F, row.names = F, col.names = T)

# process protease results


protease_pathogenic_transition <- read.delim("proteasefamily_singlelambda_k1_cafeplotter/result_summary.tsv", header = T, na.strings='', sep="\t") %>%
  filter(TaxonID %in% c(pathogenic_families_taxa)) %>%
  arrange(Pvalue) %>%
  filter(Pvalue < 0.05) %>%
  #filter(Change < 0) %>%
  left_join(orthogroup_ann, by = join_by(FamilyID == cluster_id)) %>%
  select(-protein_count, -taxon_count)
write.table(protease_pathogenic_transition, file="~/Documents/2023_02_Sporothrix/results/tables/pathogenic_transition_protease.txt", sep='\t', quote = F, row.names = F, col.names = T)

node_data <- read.delim("orthogroups_singlelambda_k1/Base_clade_results.txt", header=T) %>%
  as.data.frame() %>%
  mutate(node = str_replace(X.Taxon_ID, "(^.*<)(\\d*)(>)","\\2")) %>%
  select(-X.Taxon_ID) %>%
  relocate(node)


library("ggVennDiagram")




og_list_decrease <- list(epigloea = orthogroup_transition %>% filter(Change < 0) %>% filter(TaxonID %in% epigloea_orthogroup_taxa) %>% select(FamilyID) %>% unique() %>% pull(),
                         pathogenic = orthogroup_transition %>% filter(Change < 0) %>% filter(TaxonID %in% pathogenic_orthogroup_taxa) %>% select(FamilyID) %>% unique() %>% pull())
og_plot_decrease <- ggVennDiagram(og_list_decrease, label_alpha = 0) + 
  scale_fill_distiller(palette = "PuBu", direction = 1) +
  scale_x_continuous(expand = expansion(mult = c(0.3,0.1))) +
  labs(title = "Orthogroups\n") +
  theme(plot.title = element_text(hjust = 0.5))

cazyme_list_decrease <- list(epigloea = cazyme_transition %>% filter(Change < 0) %>% filter(TaxonID %in% epigloea_families_taxa) %>% select(FamilyID) %>% unique() %>% pull(),
                             pathogenic = cazyme_transition %>% filter(Change < 0) %>% filter(TaxonID %in% pathogenic_families_taxa) %>% select(FamilyID) %>% unique() %>% pull())
cazyme_plot_decrease <- ggVennDiagram(cazyme_list_decrease, label_alpha = 0, category.names = c("", "")) +
  scale_fill_distiller(palette = "PuBu", direction = 1) +
  scale_x_continuous(expand = expansion(mult = c(0.3,0.1))) +
  labs(title = "CAZyme\nfamilies") +
  theme(plot.title = element_text(hjust = 0.5))

protease_list_decrease <- list(epigloea = protease_transition %>% filter(Change < 0) %>% filter(TaxonID %in% epigloea_families_taxa) %>% select(FamilyID) %>% unique() %>% pull(),
                               pathogenic = protease_transition %>% filter(Change < 0) %>% filter(TaxonID %in% pathogenic_families_taxa) %>% select(FamilyID) %>% unique() %>% pull())
protease_plot_decrease <- ggVennDiagram(protease_list_decrease, label_alpha = 0, category.names = c("", "")) +
  scale_fill_distiller(palette = "PuBu", direction = 1) +
  scale_x_continuous(expand = expansion(mult = c(0.3,0.1))) +
  labs(title = "Protease\nfamilies") +
  theme(plot.title = element_text(hjust = 0.5))



og_list_increase <- list(epigloea = orthogroup_transition %>% filter(Change > 0) %>% filter(TaxonID %in% epigloea_orthogroup_taxa) %>% select(FamilyID) %>% unique() %>% pull(),
                         pathogenic = orthogroup_transition %>% filter(Change > 0) %>% filter(TaxonID %in% pathogenic_orthogroup_taxa) %>% select(FamilyID) %>% unique() %>% pull())
og_plot_increase <- ggVennDiagram(og_list_increase, label_alpha = 0) + 
  scale_fill_distiller(palette = "OrRd", direction = 1) +
  scale_x_continuous(expand = expansion(mult = c(0.3,0.1))) +
  labs(title = "Orthogroups\n") +
  theme(plot.title = element_text(hjust = 0.5))

cazyme_list_increase <- list(epigloea = cazyme_transition %>% filter(Change > 0) %>% filter(TaxonID %in% epigloea_families_taxa) %>% select(FamilyID) %>% unique() %>% pull(),
                             pathogenic = cazyme_transition %>% filter(Change > 0) %>% filter(TaxonID %in% pathogenic_families_taxa) %>% select(FamilyID) %>% unique() %>% pull())
cazyme_plot_increase <- ggVennDiagram(cazyme_list_increase, label_alpha = 0, category.names = c("", "")) +
  scale_fill_distiller(palette = "OrRd", direction = 1) +
  scale_x_continuous(expand = expansion(mult = c(0.3,0.1))) +
  labs(title = "CAZyme\nfamilies") +
  theme(plot.title = element_text(hjust = 0.5))

protease_list_increase <- list(epigloea = protease_transition %>% filter(Change > 0) %>% filter(TaxonID %in% epigloea_families_taxa) %>% select(FamilyID) %>% unique() %>% pull(),
                               pathogenic = protease_transition %>% filter(Change > 0) %>% filter(TaxonID %in% pathogenic_families_taxa) %>% select(FamilyID) %>% unique() %>% pull())
protease_plot_increase <- ggVennDiagram(protease_list_increase, label_alpha = 0, category.names = c("", "")) +
  scale_fill_distiller(palette = "OrRd", direction = 1) +
  scale_x_continuous(expand = expansion(mult = c(0.3,0.1))) +
  labs(title = "Protease\nfamilies") +
  theme(plot.title = element_text(hjust = 0.5))




pdf(file="~/Documents/2023_02_Sporothrix/results/figures/pathogenic_epigloea_overlaps.pdf", width=9, height=10)
plot_grid(og_plot_decrease, cazyme_plot_decrease, protease_plot_decrease, og_plot_increase, cazyme_plot_increase, protease_plot_increase, ncol = 3)
dev.off()


help(ggVennDiagram)

intersect(og_list_decrease[[1]],og_list_decrease[[2]])

intersect(cazyme_list_decrease[[1]],cazyme_list_decrease[[2]])

intersect(protease_list_decrease[[1]],protease_list_decrease[[2]])

intersect(og_list_increase[[1]],og_list_increase[[2]])

intersect(cazyme_list_increase[[1]],cazyme_list_increase[[2]])

intersect(protease_list_increase[[1]],protease_list_increase[[2]])


og_intersect_decrease <- intersect(og_list_decrease[[1]],og_list_decrease[[2]])
og_intersect_decrease <- as.data.frame(og_intersect_decrease, col.names = c("cluster_id"))



# not used in the end ####

######## extract the gene families that were expanded or contracted for a specific node####### used this one ####

pull_epigloea_og_contr <- function(target_gene_family) {
  
Base_change <- read.delim(paste0(target_gene_family,"_singlelambda_k1/base_change.tab", sep = ""), header = T, na.strings='', sep="\t")

Node43<- Base_change %>%
  filter(X.43. < 0) %>%
  select(FamilyID)

Node44 <- Base_change %>%
  filter(X.44. < 0) %>%
  select(FamilyID)

Node18 <- Base_change %>%
  filter(Sporothrix_epigloea_CBS57363.18. < 0) %>%
  select(FamilyID)

Node19 <- Base_change %>%
  filter(Sporothrix_epigloea_TF4163.19. < 0) %>%
  select(FamilyID)

Node20 <- Base_change %>%
  filter(Sporothrix_epigloea_CBS119000.20. < 0) %>%
  select(FamilyID)

combined <- rbind(Node43, Node44, Node18, Node19, Node20) %>% unique()

return(combined)

}

epigloea_og_contr <- pull_epigloea_og_contr(target_gene_family = "orthogroups")

pull_epigloea_og_expan <- function(target_gene_family) {
  
  Base_change <- read.delim(paste0(target_gene_family,"_singlelambda_k1/base_change.tab", sep = ""), header = T, na.strings='', sep="\t")
  
  Node43 <- Base_change %>%
    filter(X.43. > 0) %>%
    select(FamilyID)
  
  Node44 <- Base_change %>%
    filter(X.44. > 0) %>%
    select(FamilyID)
  
  Node18 <- Base_change %>%
    filter(Sporothrix_epigloea_CBS57363.18. > 0) %>%
    select(FamilyID)
  
  Node19 <- Base_change %>%
    filter(Sporothrix_epigloea_TF4163.19. > 0) %>%
    select(FamilyID)
  
  Node20 <- Base_change %>%
    filter(Sporothrix_epigloea_CBS119000.20. > 0) %>%
    select(FamilyID)
  
  combined <- rbind(Node43, Node44, Node18, Node19, Node20) %>% unique()
  
  return(combined)
  
}

epigloea_og_expan <- pull_epigloea_og_expan(target_gene_family = "orthogroups")


pull_epigloea_gfs_contr <- function(target_gene_family) {
  
  Base_change <- read.delim(paste0(target_gene_family,"_singlelambda_k1/base_change.tab", sep = ""), header = T, na.strings='', sep="\t")
  
  Node43 <- Base_change %>%
    filter(X.54. < 0) %>%
    select(FamilyID)
  
  Node44 <- Base_change %>%
    filter(X.55. < 0) %>%
    select(FamilyID)
  
  Node18 <- Base_change %>%
    filter(Sporothrix_epigloea_CBS57363.18. < 0) %>%
    select(FamilyID)
  
  Node19 <- Base_change %>%
    filter(Sporothrix_epigloea_TF4163.19. < 0) %>%
    select(FamilyID)
  
  Node20 <- Base_change %>%
    filter(Sporothrix_epigloea_CBS119000.20. < 0) %>%
    select(FamilyID)
  
  combined <- rbind(Node43, Node44, Node18, Node19, Node20) %>% unique()
  
  return(combined)
  
}

epigloea_cazyme_contr <- pull_epigloea_gfs_contr(target_gene_family = "cazymefamily")
epigloea_protease_contr <- pull_epigloea_gfs_contr(target_gene_family = "proteasefamily")

pull_epigloea_gfs_expan <- function(target_gene_family) {
  
  Base_change <- read.delim(paste0(target_gene_family,"_singlelambda_k1/base_change.tab", sep = ""), header = T, na.strings='', sep="\t")
  
  Node43 <- Base_change %>%
    filter(X.54. > 0) %>%
    select(FamilyID)
  
  Node44 <- Base_change %>%
    filter(X.55. > 0) %>%
    select(FamilyID)
  
  Node18 <- Base_change %>%
    filter(Sporothrix_epigloea_CBS57363.18. > 0) %>%
    select(FamilyID)
  
  Node19 <- Base_change %>%
    filter(Sporothrix_epigloea_TF4163.19. > 0) %>%
    select(FamilyID)
  
  Node20 <- Base_change %>%
    filter(Sporothrix_epigloea_CBS119000.20. > 0) %>%
    select(FamilyID)
  
  combined <- rbind(Node43, Node44, Node18, Node19, Node20) %>% unique()
  
  return(combined)
  
}

epigloea_cazyme_expan <- pull_epigloea_gfs_expan(target_gene_family = "cazymefamily")
epigloea_protease_expan <- pull_epigloea_gfs_expan(target_gene_family = "proteasefamily")

pull_pathogenic_og_contr <- function(target_gene_family) {
  Base_change <- read.delim(paste0(target_gene_family,"_singlelambda_k1/base_change.tab", sep = ""), header = T, na.strings='', sep="\t")
  
  Node34 <- Base_change %>%
    filter(X.34. < 0) %>%
    select(FamilyID)
  
  Node35 <- Base_change %>%
    filter(X.35. < 0) %>%
    select(FamilyID)
  
  Node36 <- Base_change %>%
    filter(X.36. < 0) %>%
    select(FamilyID)
  
  Node8 <- Base_change %>%
    filter(Sporothrix_brasiliensis_5110.8. < 0) %>%
    select(FamilyID)
  
  Node9 <- Base_change %>%
    filter(Sporothrix_schenckii_1099.9. < 0) %>%
    select(FamilyID)
  
  Node10 <- Base_change %>%
    filter(Sporothrix_globosa_CBS120340.10. < 0) %>%
    select(FamilyID)
  
  Node11 <- Base_change %>%
    filter(Sporothrix_luriei_CBS93772.11. < 0) %>%
    select(FamilyID)
  
  combined <- rbind(Node34, Node35, Node36, Node8, Node9, Node10, Node11) %>% unique()
  
  return(combined)
  
}

pathogenic_og_contr <- pull_pathogenic_og_contr(target_gene_family = "orthogroups")

pull_pathogenic_og_expan <- function(target_gene_family) {
  Base_change <- read.delim(paste0(target_gene_family,"_singlelambda_k1/base_change.tab", sep = ""), header = T, na.strings='', sep="\t")
  
  Node34 <- Base_change %>%
    filter(X.34. > 0) %>%
    select(FamilyID)
  
  Node35 <- Base_change %>%
    filter(X.35. > 0) %>%
    select(FamilyID)
  
  Node36 <- Base_change %>%
    filter(X.36. > 0) %>%
    select(FamilyID)
  
  Node8 <- Base_change %>%
    filter(Sporothrix_brasiliensis_5110.8. > 0) %>%
    select(FamilyID)
  
  Node9 <- Base_change %>%
    filter(Sporothrix_schenckii_1099.9. > 0) %>%
    select(FamilyID)
  
  Node10 <- Base_change %>%
    filter(Sporothrix_globosa_CBS120340.10. > 0) %>%
    select(FamilyID)
  
  Node11 <- Base_change %>%
    filter(Sporothrix_luriei_CBS93772.11. > 0) %>%
    select(FamilyID)
  
  combined <- rbind(Node34, Node35, Node36, Node8, Node9, Node10, Node11) %>% unique()
  
  return(combined)
  
}

pathogenic_og_expan <- pull_pathogenic_og_expan(target_gene_family = "orthogroups")

pull_pathogenic_gfs_contr <- function(target_gene_family) {
  Base_change <- read.delim(paste0(target_gene_family,"_singlelambda_k1/base_change.tab", sep = ""), header = T, na.strings='', sep="\t")
  
  Node34 <- Base_change %>%
    filter(X.45. < 0) %>%
    select(FamilyID)
  
  Node35 <- Base_change %>%
    filter(X.46. < 0) %>%
    select(FamilyID)
  
  Node36 <- Base_change %>%
    filter(X.47. < 0) %>%
    select(FamilyID)
  
  Node8 <- Base_change %>%
    filter(Sporothrix_brasiliensis_5110.8. < 0) %>%
    select(FamilyID)
  
  Node9 <- Base_change %>%
    filter(Sporothrix_schenckii_1099.9. < 0) %>%
    select(FamilyID)
  
  Node10 <- Base_change %>%
    filter(Sporothrix_globosa_CBS120340.10. < 0) %>%
    select(FamilyID)
  
  Node11 <- Base_change %>%
    filter(Sporothrix_luriei_CBS93772.11. < 0) %>%
    select(FamilyID)
  
  combined <- rbind(Node34, Node35, Node36, Node8, Node9, Node10, Node11) %>% unique()
  
  return(combined)
  
}

pathogenic_cazyme_contr <- pull_pathogenic_gfs_contr(target_gene_family = "cazymefamily")
pathogenic_protease_contr <- pull_pathogenic_gfs_contr(target_gene_family = "proteasefamily")

pull_pathogenic_gfs_expan <- function(target_gene_family) {
  Base_change <- read.delim(paste0(target_gene_family,"_singlelambda_k1/base_change.tab", sep = ""), header = T, na.strings='', sep="\t")
  
  Node34 <- Base_change %>%
    filter(X.45. > 0) %>%
    select(FamilyID)
  
  Node35 <- Base_change %>%
    filter(X.46. > 0) %>%
    select(FamilyID)
  
  Node36 <- Base_change %>%
    filter(X.47. > 0) %>%
    select(FamilyID)
  
  Node8 <- Base_change %>%
    filter(Sporothrix_brasiliensis_5110.8. > 0) %>%
    select(FamilyID)
  
  Node9 <- Base_change %>%
    filter(Sporothrix_schenckii_1099.9. > 0) %>%
    select(FamilyID)
  
  Node10 <- Base_change %>%
    filter(Sporothrix_globosa_CBS120340.10. > 0) %>%
    select(FamilyID)
  
  Node11 <- Base_change %>%
    filter(Sporothrix_luriei_CBS93772.11. > 0) %>%
    select(FamilyID)
  
  combined <- rbind(Node34, Node35, Node36, Node8, Node9, Node10, Node11) %>% unique()
  
  return(combined)
  
}

pathogenic_cazyme_expan <- pull_pathogenic_gfs_expan(target_gene_family = "cazymefamily")
pathogenic_protease_expan <- pull_pathogenic_gfs_expan(target_gene_family = "proteasefamily")

library("ggVennDiagram")



og_list_contr <- list(epigloea = epigloea_og_contr %>% pull(),
                         pathogenic = pathogenic_og_contr %>% pull())
og_plot_contr <- ggVennDiagram(og_list_contr, label_alpha = 0, category.names = c("epigloea", "pathogenic")) + 
  scale_fill_distiller(palette = "PuBu", direction = 1) +
  scale_x_continuous(expand = expansion(mult = c(0.3,0.1))) +
  labs(title = "Orthogroups\n") +
  theme(plot.title = element_text(hjust = 0.5))

cazyme_list_contr <- list(epigloea = epigloea_cazyme_contr %>% pull(),
                      pathogenic = pathogenic_cazyme_contr %>% pull())
cazyme_plot_contr <- ggVennDiagram(cazyme_list_contr, label_alpha = 0, category.names = c("", "")) + 
  scale_fill_distiller(palette = "PuBu", direction = 1) +
  scale_x_continuous(expand = expansion(mult = c(0.3,0.1))) +
  labs(title = "CAZyme\nfamilies") +
  theme(plot.title = element_text(hjust = 0.5))

protease_list_contr <- list(epigloea = epigloea_protease_contr %>% pull(),
                          pathogenic = pathogenic_protease_contr %>% pull())
protease_plot_contr <- ggVennDiagram(protease_list_contr, label_alpha = 0, category.names = c("", "")) + 
  scale_fill_distiller(palette = "PuBu", direction = 1) +
  scale_x_continuous(expand = expansion(mult = c(0.3,0.1))) +
  labs(title = "Protease\nfamilies") +
  theme(plot.title = element_text(hjust = 0.5))

og_list_expan <- list(epigloea = epigloea_og_expan %>% pull(),
                      pathogenic = pathogenic_og_expan %>% pull())
og_plot_expan <- ggVennDiagram(og_list_expan, label_alpha = 0, category.names = c("epigloea", "pathogenic")) + 
  scale_fill_distiller(palette = "OrRd", direction = 1) +
  scale_x_continuous(expand = expansion(mult = c(0.3,0.1))) +
  #labs(title = "Orthogroups\n") +
  theme(plot.title = element_text(hjust = 0.5))

cazyme_list_expan <- list(epigloea = epigloea_cazyme_expan %>% pull(),
                          pathogenic = pathogenic_cazyme_expan %>% pull())
cazyme_plot_expan <- ggVennDiagram(cazyme_list_expan, label_alpha = 0, category.names = c("", "")) + 
  scale_fill_distiller(palette = "OrRd", direction = 1) +
  scale_x_continuous(expand = expansion(mult = c(0.3,0.1))) +
  #labs(title = "CAzyme\nfamilies") +
  theme(plot.title = element_text(hjust = 0.5))

protease_list_expan <- list(epigloea = epigloea_protease_expan %>% pull(),
                            pathogenic = pathogenic_protease_expan %>% pull())
protease_plot_expan <- ggVennDiagram(protease_list_expan, label_alpha = 0, category.names = c("", "")) + 
  scale_fill_distiller(palette = "OrRd", direction = 1) +
  scale_x_continuous(expand = expansion(mult = c(0.3,0.1))) +
  #labs(title = "Protease\nfamilies") +
  theme(plot.title = element_text(hjust = 0.5))

pdf(file="~/Documents/2023_02_Sporothrix/results/figures/pathogenic_epigloea_overlaps.pdf", width=9, height=10)
plot_grid(og_plot_contr, cazyme_plot_contr, protease_plot_contr, og_plot_expan, cazyme_plot_expan, protease_plot_expan, ncol = 3)
dev.off()

#################################################################

#extract intersections


intersect(og_list_contr[[1]],og_list_contr[[2]])
intersect(cazyme_list_contr[[1]],cazyme_list_contr[[2]])
intersect(protease_list_contr[[1]],protease_list_contr[[2]])

intersect(og_list_expan[[1]],og_list_expan[[2]])
intersect(cazyme_list_expan[[1]],cazyme_list_expan[[2]])
intersect(protease_list_expan[[1]],protease_list_expan[[2]])

setdiff(og_list_contr[[1]],og_list_contr[[2]])
setdiff(cazyme_list_contr[[1]],cazyme_list_contr[[2]])
setdiff(protease_list_contr[[1]],protease_list_contr[[2]])

setdiff(og_list_expan[[1]],og_list_expan[[2]])
setdiff(cazyme_list_expan[[1]],cazyme_list_expan[[2]])
setdiff(protease_list_expan[[1]],protease_list_expan[[2]])

setdiff(og_list_contr[[2]],og_list_contr[[1]])
setdiff(cazyme_list_contr[[2]],cazyme_list_contr[[1]])
setdiff(protease_list_contr[[2]],protease_list_contr[[1]])

setdiff(og_list_expan[[2]],og_list_expan[[1]])
setdiff(cazyme_list_expan[[2]],cazyme_list_expan[[1]])
setdiff(protease_list_expan[[2]],protease_list_expan[[1]])




draw_intersect_wordcloud <- function(target_set) {

  intersect <- intersect(target_set[[1]],target_set[[2]])
  intersect <- as.data.frame(intersect)
  colnames(intersect) = c("cluster_id")
  intersect <- intersect %>%
  left_join(orthogroup_ann) %>%
  select(IPR_domain_description) %>%
  separate_longer_delim(IPR_domain_description, delim = ";") %>%
  filter(IPR_domain_description != "None") %>%
  group_by(IPR_domain_description) %>%
  count() %>%
  ungroup()
  colnames(intersect) = c("word", "freq")
  
  set.seed(1234)
  wordcloud <- wordcloud(words = intersect$word,
          freq = intersect$freq,
          scale=c(2,.5),
          min.freq = 4,
          max.words=100,
          random.order=FALSE,
          rot.per=0, #	proportion words with 90 degree rotation
          colors=brewer.pal(8, "Dark2"))
  return(wordcloud)
}
draw_diffepigloea_wordcloud <- function(target_set) {
  
  diff <- setdiff(target_set[[1]],target_set[[2]])
  diff <- as.data.frame(diff)
  colnames(diff) = c("cluster_id")
  diff <- diff %>%
    left_join(orthogroup_ann) %>%
    select(IPR_domain_description) %>%
    separate_longer_delim(IPR_domain_description, delim = ";") %>%
    filter(IPR_domain_description != "None") %>%
    group_by(IPR_domain_description) %>%
    count() %>%
    ungroup()
  colnames(diff) = c("word", "freq")
  
  set.seed(1234)
  wordcloud <- wordcloud(words = diff$word,
                         freq = diff$freq,
                         scale=c(2,.5),
                         min.freq = 4,
                         max.words=100,
                         random.order=FALSE,
                         rot.per=0, #	proportion words with 90 degree rotation
                         colors=brewer.pal(8, "Dark2"))
  return(wordcloud)
}
draw_diffpathogenic_wordcloud <- function(target_set) {
  
  diff <- setdiff(target_set[[2]],target_set[[1]])
  diff <- as.data.frame(diff)
  colnames(diff) = c("cluster_id")
  diff <- diff %>%
    left_join(orthogroup_ann) %>%
    select(IPR_domain_description) %>%
    separate_longer_delim(IPR_domain_description, delim = ";") %>%
    filter(IPR_domain_description != "None") %>%
    group_by(IPR_domain_description) %>%
    count() %>%
    ungroup()
  colnames(diff) = c("word", "freq")
  
  set.seed(1234)
  wordcloud <- wordcloud(words = diff$word,
                         freq = diff$freq,
                         scale=c(2,.5),
                         min.freq = 4,
                         max.words=100,
                         random.order=FALSE,
                         rot.per=0, #	proportion words with 90 degree rotation
                         colors=brewer.pal(8, "Dark2"))
  return(wordcloud)
}
draw_intersect_wordcloud(og_list_contr)
draw_intersect_wordcloud(og_list_expan)
draw_diffepigloea_wordcloud(og_list_contr)
draw_diffepigloea_wordcloud(og_list_expan)
draw_diffpathogenic_wordcloud(og_list_contr)
draw_diffpathogenic_wordcloud(og_list_expan)


pdf(file="~/Documents/2023_02_Sporothrix/results/figures/wordcloud_og_shared_down.pdf", width=10, height=10)
draw_intersect_wordcloud(og_list_contr)
dev.off()

pdf(file="~/Documents/2023_02_Sporothrix/results/figures/wordcloud_og_epigloea_down.pdf", width=10, height=10)
draw_diffepigloea_wordcloud(og_list_contr)
dev.off()

pdf(file="~/Documents/2023_02_Sporothrix/results/figures/wordcloud_og_pathogen_down.pdf", width=10, height=10)
draw_diffpathogenic_wordcloud(og_list_contr)
dev.off()

pdf(file="~/Documents/2023_02_Sporothrix/results/figures/wordcloud_og_shared_up.pdf", width=10, height=10)
draw_intersect_wordcloud(og_list_expan)
dev.off()

pdf(file="~/Documents/2023_02_Sporothrix/results/figures/wordcloud_og_epigloea_up.pdf", width=10, height=10)
draw_diffepigloea_wordcloud(og_list_expan)
dev.off()

pdf(file="~/Documents/2023_02_Sporothrix/results/figures/wordcloud_og_pathogen_up.pdf", width=10, height=10)
draw_diffpathogenic_wordcloud(og_list_expan)
dev.off()


epigloea_og_expan_ann <- epigloea_og_expan
colnames(epigloea_og_expan_ann) = c("cluster_id")
epigloea_og_expan_ann <- epigloea_og_expan_ann %>%
  left_join(orthogroup_ann) %>%
  select(IPR_domain_description) %>%
  separate_longer_delim(IPR_domain_description, delim = ";") %>%
  filter(IPR_domain_description != "None") %>%
  group_by(IPR_domain_description) %>%
  count() %>%
  ungroup()
colnames(epigloea_og_expan_ann) = c("word", "freq")



#####################################################################
# extract orthogroups that increased or decreased at the transition branch to S. epigloea
  
Base_change <- read.delim("orthogroups_singlelambda_k1/base_change.tab", header = T, na.strings='', sep="\t")
Node43og <- Base_change %>%
    #filter(X.43. < 0) %>%
  filter(X.43. != 0) %>%
  select(FamilyID, X.43.) %>%
  arrange(X.43.)
colnames(Node43og) = c("cluster_id", "change_node43")
Node43og <- Node43og %>%
  left_join(orthogroup_ann)  %>%
  select(-protein_count,-taxon_count)

Node43down <- Node43og %>% filter(change_node43 < 0)
write.table(Node43down, file="~/Documents/2023_02_Sporothrix/results/tables/epigloea_orthogroup_contraction_annotations.txt", sep='\t', quote = F, row.names = F, col.names = T)
Node43up <- Node43og %>% filter(change_node43 > 0) %>% arrange(desc(change_node43))
write.table(Node43up, file="~/Documents/2023_02_Sporothrix/results/tables/epigloea_orthogroup_expansion_annotations.txt", sep='\t', quote = F, row.names = F, col.names = T)


Node43downwordcloud <- Node43down %>%
  select(IPR_domain_description) %>%
  separate_longer_delim(IPR_domain_description, delim = ";") %>%
  filter(IPR_domain_description != "None") %>%
  group_by(IPR_domain_description) %>%
  count() %>%
  ungroup()
colnames(Node43downwordcloud) = c("word", "freq")

pdf(file="~/Documents/2023_02_Sporothrix/results/figures/wordcloud_epigloea_contr.pdf", width=10, height=10)
set.seed(1234)
wordcloud <- wordcloud(words = Node43downwordcloud$word,
                       freq = Node43downwordcloud$freq,
                       scale=c(2,.5),
                       min.freq = 4,
                       max.words=100,
                       random.order=FALSE,
                       rot.per=0, #	proportion words with 90 degree rotation
                       colors=brewer.pal(8, "Dark2"))
dev.off()

# extract and count the IPR annotations

#descriptions
Node43downcount <- Node43down %>%
  select(IPR_domain_description) %>%
  separate_longer_delim(IPR_domain_description, delim = ";") %>%
  #filter(IPR_domain_ids != "None") %>%
  group_by(IPR_domain_description) %>%
  count() %>%
  ungroup() %>%
  arrange(desc(n))

write.table(Node43downcount, file="~/Documents/2023_02_Sporothrix/results/tables/epigloea_orthogroup_contraction_annotations_count.txt", sep='\t', quote = F, row.names = F, col.names = T)

Node43upcount <- Node43up %>%
  select(IPR_domain_description) %>%
  separate_longer_delim(IPR_domain_description, delim = ";") %>%
  #filter(IPR_domain_ids != "None") %>%
  group_by(IPR_domain_description) %>%
  count() %>%
  ungroup() %>%
  arrange(desc(n))

write.table(Node43upcount, file="~/Documents/2023_02_Sporothrix/results/tables/epigloea_orthogroup_expansion_annotations_count.txt", sep='\t', quote = F, row.names = F, col.names = T)

####################################

# extract CAZymes that increased or decreased at the transition branch to S. epigloea

Base_change <- read.delim("cazymefamily_singlelambda_k1/base_change.tab", header = T, na.strings='', sep="\t")
Node54cazyme <- Base_change %>%
  #filter(X.54. < 0) %>%
  filter(X.54. != 0) %>%
  select(FamilyID, X.54.) %>% arrange(X.54.)
colnames(Node54cazyme) = c("FamilyID", "change_node54")
write.table(Node54cazyme, file="~/Documents/2023_02_Sporothrix/results/tables/epigloea_cazyme_node54change.txt", sep='\t', quote = F, row.names = F, col.names = T)


Node54cazymewordcloud <- Node54cazyme
Node54cazymewordcloud$freq <- abs(Node54cazymewordcloud$change_node54)
Node54cazymewordcloud <- Node54cazymewordcloud %>% select(-change_node54)
colnames(Node54cazymewordcloud) = c("word", "freq")

pdf(file="~/Documents/2023_02_Sporothrix/results/figures/wordcloud_epigloea_cazyme_contr.pdf", width=10, height=10)
set.seed(1234)
wordcloud <- wordcloud(words = Node54cazymewordcloud$word,
                       freq = Node54cazymewordcloud$freq,
                       scale=c(2,.5),
                       min.freq = 1,
                       max.words=100,
                       random.order=FALSE,
                       rot.per=0, #	proportion words with 90 degree rotation
                       colors=brewer.pal(8, "Dark2"))
dev.off()


####################################

# extract proteases that increased or decreased at the transition branch to S. epigloea

Base_change <- read.delim("proteasefamily_singlelambda_k1/base_change.tab", header = T, na.strings='', sep="\t")
Node54protease <- Base_change %>%
  #filter(X.54. < 0) %>%
  filter(X.54. != 0) %>%
  select(FamilyID, X.54.) %>% arrange(X.54.)
colnames(Node54protease) = c("FamilyID", "change_node54")
write.table(Node54protease, file="~/Documents/2023_02_Sporothrix/results/tables/epigloea_protease_node54change.txt", sep='\t', quote = F, row.names = F, col.names = T)


Node54proteasewordcloud <- Node54protease %>% filter(change_node54 < 0)
Node54proteasewordcloud$freq <- abs(Node54proteasewordcloud$change_node54)
Node54proteasewordcloud <- Node54proteasewordcloud %>% select(-change_node54)
colnames(Node54proteasewordcloud) = c("word", "freq")

pdf(file="~/Documents/2023_02_Sporothrix/results/figures/wordcloud_epigloea_protease_contr.pdf", width=10, height=10)
set.seed(1234)
wordcloud <- wordcloud(words = Node54proteasewordcloud$word,
                       freq = Node54proteasewordcloud$freq,
                       scale=c(2,.5),
                       min.freq = 1,
                       max.words=100,
                       random.order=FALSE,
                       rot.per=0, #	proportion words with 90 degree rotation
                       colors=brewer.pal(8, "Dark2"))
dev.off()

