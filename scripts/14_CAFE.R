
# make an ultrametric tree for CAFE analysis

setwd("~/Documents/2023_02_Sporothrix/14_CAFE/")
library(ape)
library(phytools)

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
write.tree(ML_tree, file="~/Documents/2023_02_Sporothrix/09_orthology/OrthoFinder/Results_Aug03/Species_Tree/iqtree2/sporothrix_placement.rerooted.ultrametric.contree")

# format protease family results from funannotate for CAFE5

protease_family <- read.delim("~/Documents/2023_02_Sporothrix/10_annotations/funannotate/funannotate_compare_20230809/merops/MEROPS.all.results.csv", header = T, na.strings='', sep=",") %>%
  mutate(Desc = "(null)") %>%
  relocate(Desc)
colnames(protease_family)[2] <- "Family ID"
write.table(protease_family, file="~/Documents/2023_02_Sporothrix/14_CAFE/protease_families_funannotate.txt", sep='\t', quote = F, row.names = F, col.names = T)

# format protease class results from funannotate for CAFE5
# in the end I did not use this

protease_class <- read.delim("~/Documents/2023_02_Sporothrix/10_annotations/funannotate/funannotate_compare_20230809/merops/MEROPS.summary.results.csv", header = T, na.strings='', sep=",") %>%
  mutate(Desc = "(null)") %>%
  relocate(Desc)
colnames(protease_class)[2] <- "Family ID"
write.table(protease_class, file="~/Documents/2023_02_Sporothrix/14_CAFE/protease_classes_funannotate.txt", sep='\t', quote = F, row.names = F, col.names = T)

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
write.table(secmets, file="~/Documents/2023_02_Sporothrix/14_CAFE/secmets.txt", sep='\t', quote = F, row.names = F, col.names = T)

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
write.tree(ML_tree_orthog, file="~/Documents/2023_02_Sporothrix/09_orthology/OrthoFinder/Results_Aug03/Species_Tree/iqtree2/sporothrix_placement.rerooted.ultrametric.coresporothrix.contree")

orthogroups <- read.delim("~/Documents/2023_02_Sporothrix/09_orthology/OrthoFinder/Results_Aug07/Orthogroups/Orthogroups.GeneCount.tsv", header = T, na.strings='', sep="\t") %>%
  mutate(Desc = "(null)") %>%
  relocate(Desc) %>%
  select(-Total)
colnames(orthogroups)[2] <- "Family ID"
write.table(orthogroups, file="~/Documents/2023_02_Sporothrix/14_CAFE/orthogroups.txt", sep='\t', quote = F, row.names = F, col.names = T)

# visualize CAFE results

install.packages("ggtree")

library(tidyverse)
library(ggtree)
library(treeio)
library(tidytree)



tree <- read.tree("~/Documents/2023_02_Sporothrix/09_orthology/OrthoFinder/Results_Aug03/Species_Tree/iqtree2/sporothrix_placement.rerooted.ultrametric.contree")

#draw tree with node numbers
#ggtree(tree) + geom_text(aes(label=node), hjust=-.3) + geom_tiplab(hjust=-0.3)

# convert tree into tibble tree
tbl_tree <- as_tibble(tree)
as.treedata(tbl_tree)
ggtree(tbl_tree)
str(tbl_tree)

# extract number of families from CAFE analysis
node_data <- read.delim("cazymefamily_singlelambda_k1/Base_clade_results.txt", header=T) %>%
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

#child(tree_data, 51)


ggtree(tree_data) +
  geom_tiplab(offset=.05) +
  theme_tree2() +
  ggplot2::xlim(0, 2) +
  geom_point(aes(size = Decrease)) +
  scale_size(range = c(0, 11))





# use the above to create plotting functions

plot_CAFE_decrease <- function(target_tree, target_data, target_group){

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

p <- ggtree(tree_data) +
  geom_tiplab(offset=.05, size = 3) +
  theme_tree() +
  ggplot2::xlim(0, 2) +
  vexpand(.05, direction = -1) +
  geom_point(aes(size = Decrease), colour="skyblue4", alpha=0.5) +
  scale_size(range = c(0, 11)) +
  ggtitle(paste("Decrease - ",target_group,sep = "")) +
  theme(plot.title = element_text(hjust = 0.5))

return(p)

}

plot_CAFE_increase <- function(target_tree, target_data, target_group){
  
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
  
  p <- ggtree(tree_data) +
    geom_tiplab(offset=.05, size =3) +
    theme_tree() +
    ggplot2::xlim(0, 2) +
    vexpand(.05, direction = -1) +
    geom_point(aes(size = Increase), colour="khaki3", alpha=0.5) +
    scale_size(range = c(0, 11)) +
    ggtitle(paste("Increase - ",target_group,sep = "")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
  
}

cazyme_decrease <- plot_CAFE_decrease("~/Documents/2023_02_Sporothrix/09_orthology/OrthoFinder/Results_Aug03/Species_Tree/iqtree2/sporothrix_placement.rerooted.ultrametric.contree",
                                      "cazymefamily_singlelambda_k1/Base_clade_results.txt", target_group = "CAZyme families")

protease_decrease <- plot_CAFE_decrease("~/Documents/2023_02_Sporothrix/09_orthology/OrthoFinder/Results_Aug03/Species_Tree/iqtree2/sporothrix_placement.rerooted.ultrametric.contree",
                                        "proteasefamily_singlelambda_k1/Base_clade_results.txt", target_group = "Protease families")

secmets_decrease <- plot_CAFE_decrease("~/Documents/2023_02_Sporothrix/09_orthology/OrthoFinder/Results_Aug03/Species_Tree/iqtree2/sporothrix_placement.rerooted.ultrametric.contree",
                                        "secmets_singlelambda_k6/Gamma_clade_results.txt", target_group = "BGCs")

orthogroups_decrease <- plot_CAFE_decrease("~/Documents/2023_02_Sporothrix/09_orthology/OrthoFinder/Results_Aug03/Species_Tree/iqtree2/sporothrix_placement.rerooted.ultrametric.coresporothrix.contree",
                                       "orthogroups_singlelambda_k1/Base_clade_results.txt", target_group = "Orthogroups")

cazyme_increase <- plot_CAFE_increase("~/Documents/2023_02_Sporothrix/09_orthology/OrthoFinder/Results_Aug03/Species_Tree/iqtree2/sporothrix_placement.rerooted.ultrametric.contree",
                                      "cazymefamily_singlelambda_k1/Base_clade_results.txt", target_group = "CAZyme families")

protease_increase <- plot_CAFE_increase("~/Documents/2023_02_Sporothrix/09_orthology/OrthoFinder/Results_Aug03/Species_Tree/iqtree2/sporothrix_placement.rerooted.ultrametric.contree",
                                        "proteasefamily_singlelambda_k1/Base_clade_results.txt", target_group = "Protease families")

secmets_increase <- plot_CAFE_increase("~/Documents/2023_02_Sporothrix/09_orthology/OrthoFinder/Results_Aug03/Species_Tree/iqtree2/sporothrix_placement.rerooted.ultrametric.contree",
                                       "secmets_singlelambda_k6/Gamma_clade_results.txt", target_group = "BGCs")

orthogroups_increase <- plot_CAFE_increase("~/Documents/2023_02_Sporothrix/09_orthology/OrthoFinder/Results_Aug03/Species_Tree/iqtree2/sporothrix_placement.rerooted.ultrametric.coresporothrix.contree",
                                           "orthogroups_singlelambda_k1/Base_clade_results.txt", target_group = "Orthogroups")


library(cowplot)

pdf(file="~/Documents/2023_02_Sporothrix/results/figures/families_expansion_contraction.pdf", width=12, height=10)
plot_grid(cazyme_decrease, protease_decrease, cazyme_increase, protease_increase, ncol = 2)
dev.off()


pdf(file="~/Documents/2023_02_Sporothrix/results/figures/orthogroups_expansion_contraction.pdf", width=10, height=5)
plot_grid(orthogroups_decrease, orthogroups_increase, ncol = 2)
dev.off()



pdf(file="~/Documents/2023_02_Sporothrix/results/figures/all_expansion_contraction.pdf", width=25, height=10)
plot_grid(orthogroups_decrease, cazyme_decrease, protease_decrease, secmets_decrease, orthogroups_increase, cazyme_increase, protease_increase, secmets_increase, ncol = 4)
dev.off()




# process orthogroup results

# load kinfin results
GO_rep <- read.delim("../09_orthology/OrthoFinder/KinFin/GO.cluster_functional_annotation.all.p30.x2.tsv", sep = "\t", stringsAsFactors = F, col.names = c("cluster_id", "protein_count", "taxon_count", "GO_domain_ids", "GO_domain_description"))
IPR_rep <- read.delim("../09_orthology/OrthoFinder/KinFin/IPR.cluster_functional_annotation.all.p30.x2.tsv", sep = "\t", stringsAsFactors = F, col.names = c("cluster_id", "protein_count", "taxon_count", "IPR_domain_ids", "IPR_domain_description"))
Pfam_rep <- read.delim("../09_orthology/OrthoFinder/KinFin/Pfam.cluster_functional_annotation.all.p30.x2.tsv", sep = "\t", stringsAsFactors = F, col.names = c("cluster_id", "protein_count", "taxon_count", "Pfam_domain_ids", "Pfam_domain_description"))

#combine the 3 annotation methods
orthogroup_ann <- left_join(GO_rep, IPR_rep)
orthogroup_ann <- left_join(orthogroup_ann, Pfam_rep)

orthogroup_decrease_epigloea_transition <- read.delim("orthogroups_singlelambda_k1_cafeplotter/result_summary.tsv", header = T, na.strings='', sep="\t") %>%
  filter(TaxonID == 43|TaxonID == 44|TaxonID == 18|TaxonID == 19|TaxonID == 20) %>%
  arrange(Pvalue) %>%
  filter(Pvalue < 0.05) %>%
  filter(Change < 0) %>%
  left_join(orthogroup_ann, by = join_by(FamilyID == cluster_id)) %>%
  select(-protein_count, -taxon_count)
write.table(orthogroup_decrease_epigloea_transition, file="~/Documents/2023_02_Sporothrix/results/tables/orthogroup_decrease_epigloea_transition.txt", sep='\t', quote = F, row.names = F, col.names = T)


orthogroup_increase_epigloea_transition <- read.delim("orthogroups_singlelambda_k1_cafeplotter/result_summary.tsv", header = T, na.strings='', sep="\t") %>%
  filter(TaxonID == 43|TaxonID == 44|TaxonID == 18|TaxonID == 19|TaxonID == 20) %>%
  arrange(Pvalue) %>%
  filter(Pvalue < 0.05) %>%
  filter(Change > 0) %>%
  left_join(orthogroup_ann, by = join_by(FamilyID == cluster_id)) %>%
  select(-protein_count, -taxon_count)
write.table(orthogroup_increase_epigloea_transition, file="~/Documents/2023_02_Sporothrix/results/tables/orthogroup_increase_epigloea_transition.txt", sep='\t', quote = F, row.names = F, col.names = T)

# process cazyme results


cazyme_decrease_epigloea_transition <- read.delim("cazymefamily_singlelambda_k1_cafeplotter/result_summary.tsv", header = T, na.strings='', sep="\t") %>%
  filter(TaxonID == 54|TaxonID == 55|TaxonID == 18|TaxonID == 19|TaxonID == 20) %>%
  arrange(Pvalue) %>%
  filter(Pvalue < 0.05) %>%
  filter(Change < 0)
write.table(cazyme_decrease_epigloea_transition, file="~/Documents/2023_02_Sporothrix/results/tables/cazyme_decrease_epigloea_transition.txt", sep='\t', quote = F, row.names = F, col.names = T)


cazyme_increase_epigloea_transition <- read.delim("cazymefamily_singlelambda_k1_cafeplotter/result_summary.tsv", header = T, na.strings='', sep="\t") %>%
  filter(TaxonID == 54|TaxonID == 55|TaxonID == 18|TaxonID == 19|TaxonID == 20) %>%
  arrange(Pvalue) %>%
  filter(Pvalue < 0.05) %>%
  filter(Change > 0)
write.table(cazyme_increase_epigloea_transition, file="~/Documents/2023_02_Sporothrix/results/tables/orthogroup_increase_epigloea_transition.txt", sep='\t', quote = F, row.names = F, col.names = T)


# process protease results


protease_decrease_epigloea_transition <- read.delim("proteasefamily_singlelambda_k1_cafeplotter/result_summary.tsv", header = T, na.strings='', sep="\t") %>%
  filter(TaxonID == 54|TaxonID == 55|TaxonID == 18|TaxonID == 19|TaxonID == 20) %>%
  arrange(Pvalue) %>%
  filter(Pvalue < 0.05) %>%
  filter(Change < 0)
write.table(protease_decrease_epigloea_transition, file="~/Documents/2023_02_Sporothrix/results/tables/protease_decrease_epigloea_transition.txt", sep='\t', quote = F, row.names = F, col.names = T)


protease_increase_epigloea_transition <- read.delim("proteasefamily_singlelambda_k1_cafeplotter/result_summary.tsv", header = T, na.strings='', sep="\t") %>%
  filter(TaxonID == 54|TaxonID == 55|TaxonID == 18|TaxonID == 19|TaxonID == 20) %>%
  arrange(Pvalue) %>%
  filter(Pvalue < 0.05) %>%
  filter(Change > 0)
write.table(protease_increase_epigloea_transition, file="~/Documents/2023_02_Sporothrix/results/tables/protease_increase_epigloea_transition.txt", sep='\t', quote = F, row.names = F, col.names = T)

# process secmets results



