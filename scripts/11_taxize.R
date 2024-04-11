#install.packages("taxize")
library(taxize)
library(tidyverse)
#########################

setwd("~/Documents/2023_02_Sporothrix/11_LPMO_trees/AA14/Fungi_March2024/")
headers_blast <- readLines("headers_blast.txt")

# make a tibble with tips, accesson number, and genus as extracted from the headers_blast
blast_table <- tibble(
  tip = str_extract(headers_blast, ".*"),
  accession = str_extract(headers_blast, "(^.*?\\.\\d)\\_", group=1),
  genus = str_extract(headers_blast, "\\_([A-Z].[a-z]*)\\_", group=1)
)

genome_table <- read_delim("headers_genome.txt", delim = "\t", col_names = TRUE)

# combine blast results and genome annotation results and get distinct genera
query <- blast_table %>% bind_rows(genome_table) %>% select(genus) %>% distinct()

#usethis::edit_r_environ()
#ENTREZ_KEY='d95cd98b3bf9d147137df71e5659fd67c108'
#help(getkey)
#help(use_entrez)
#taxize::use_entrez()

#run classification
list_out <- classification(query$genus, db = "ncbi")

# a tibble of the results
t1 <- rbind(list_out) %>% 
  dplyr::select(name, rank, query) %>% 
  filter(str_detect(rank, '^kingdom|^phylum|^class|^order|^family')) %>% 
  distinct(query, rank, .keep_all = TRUE) %>% 
  pivot_wider(names_from = rank,
              values_from = name) %>% 
  replace_na(list(kingdom = "unknown", phylum = "unknown", class = "unknown", order = "unknown", family = "unknown")) %>% 
  #filter(kingdom != "Fungi") %>%
  dplyr::select(kingdom, phylum, class, order, family, query) %>%
  dplyr::rename(genus = query)


# combine with the tips
t2 <- left_join(blast_table, t1, by = "genus") %>%
  arrange(kingdom, phylum, class, order, family)
  
#write.table(t2, "tree_taxonomy.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
###################
t2 <- read_delim("itol/tree_taxonomy.txt", del = "\t", col_names = TRUE)

# add colours for itol

library(RColorBrewer)
library(circlize)

col_fun = colorRamp2(c(1, length(unique(t2$kingdom))), c("gray23", "gray90"))
colours_kingdom <- col_fun(seq(1, length(unique(t2$kingdom))))
#colours_kingdom <- brewer.pal(n = length(unique(t2$kingdom)), "Dark2")
#names(colours_kingdom) = unique(t2$kingdom)
colours_kingdom <- as.data.frame(colours_kingdom)
colours_kingdom$kingdom = unique(t2$kingdom)

legend_kingdom <- colours_kingdom %>% mutate(LEGEND_SHAPES = 1, .before=colours_kingdom)
colnames(legend_kingdom) = c("LEGEND_SHAPES", "LEGEND_COLORS", "LEGEND_LABELS")
legend_kingdom = t(legend_kingdom)
write.table(legend_kingdom, "itol/kingdom_legend.txt", sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)

t3 <- left_join(t2, colours_kingdom, by = "kingdom")

fungi <- t3 %>% filter(kingdom == "Fungi")
#col_fun = colorRamp2(c(1, length(unique(fungi$phylum))), c("gray23", "#88A0DC"))
#colours_phylum <- col_fun(seq(1, length(unique(fungi$phylum))))
colours_phylum <- rev(brewer.pal(n = length(unique(fungi$phylum)), "Blues"))
colours_phylum

colours_phylum <- as.data.frame(colours_phylum)
colours_phylum$phylum = unique(fungi$phylum)

legend_phylum <- colours_phylum %>% mutate(LEGEND_SHAPES = 1, .before=colours_phylum)
colnames(legend_phylum) = c("LEGEND_SHAPES", "LEGEND_COLORS", "LEGEND_LABELS")
legend_phylum = t(legend_phylum)
write.table(legend_phylum, "itol/phylum_legend.txt", sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)

t3 <- left_join(t3, colours_phylum, by = "phylum")

# pull out Sordariomycetes and set colours for orders
sordario <- t3 %>% filter(class == "Sordariomycetes")
#colours_sordario <- brewer.pal(n = length(unique(sordario$order)), "Paired")
colours_sordario <- rev(brewer.pal(n = length(unique(sordario$order)), "Greens"))
colours_sordario <- as.data.frame(colours_sordario)
colours_sordario$order = unique(sordario$order)

legend_sordario <- colours_sordario %>% mutate(LEGEND_SHAPES = 1, .before=colours_sordario)
colnames(legend_sordario) = c("LEGEND_SHAPES", "LEGEND_COLORS", "LEGEND_LABELS")
legend_sordario = t(legend_sordario)
write.table(legend_sordario, "itol/sordario_legend.txt", sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)

t3 <- full_join(t3, colours_sordario, by = "order")

# pull out Hypocreales and set colours for orders
hypocr <- t3 %>% filter(order == "Hypocreales")
#unique(hypocr$family)
library(MetBrewer)
#colours_hypocr = met.brewer(name="Archambault", n=1+length(unique(hypocr$family)))
#colours_hypocr <- colours_hypocr[2:9]
#colours_hypocr
family <- c("Bionectriaceae", "Clavicipitaceae", "Cordycipitaceae", "Hypocreaceae", "Nectriaceae", "Ophiocordycipitaceae", "Stachybotryaceae")
colours_hypocr <- c("#cab2d6","#5A326A", "#c1cdc1","#ff7f00","#ff3030","#8b4513", "#FFFF00" )
colours_hypocr <- data.frame(colours_hypocr, family)


#colours_hypocr <- brewer.pal(n = length(unique(hypocr$family)), "Paired")
#colours_hypocr <- as.data.frame(colours_hypocr)
#colours_hypocr$family = unique(hypocr$family)
#colours_hypocr <- colours_hypocr %>% filter(family != "unknown")

legend_hypocr <- colours_hypocr %>% mutate(LEGEND_SHAPES = 1, .before=colours_hypocr)
colnames(legend_hypocr) = c("LEGEND_SHAPES", "LEGEND_COLORS", "LEGEND_LABELS")
legend_hypocr = t(legend_hypocr)
write.table(legend_hypocr, "itol/hypocr_legend.txt", sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)


t3 <- full_join(t3, colours_hypocr, by = "family")

write.table(t3, "itol/taxonomy_colours.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# extract kingdom colours

kingdom_ann <- t3 %>% select(tip, colours_kingdom)
write.table(kingdom_ann, "itol/kingdom_colours.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# extract phylum colours

phylum_ann <- t3 %>% select(tip, colours_phylum) %>%
  drop_na(colours_phylum)
write.table(phylum_ann, "itol/phylum_colours.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# extract Sordario colours

sordario_ann <- t3 %>% select(tip, colours_sordario) %>%
  drop_na(colours_sordario)
write.table(sordario_ann, "itol/sordario_colours.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


# extract Hypocreales colours

hypocreales_ann <- t3 %>% select(tip, colours_hypocr) %>%
  drop_na(colours_hypocr)
write.table(hypocreales_ann, "itol/hypocreales_colours.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


