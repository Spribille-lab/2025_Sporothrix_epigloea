rlibrary(tidyverse)
library(patchwork)
library(ggrepel)
library(data.table)

# this is how I made the big dataset

setwd("~/Documents/2022_10_phenotypes/02_data_visualization/")

options(ggrepel.max.overlaps = Inf) # allows for infinite overlaps of labels

# import plate information
PM1sources <- read.delim("PM1sources.txt", stringsAsFactors = F, sep="\t") %>%
  select(array, position, carbon_source)
PM2sources <- read.delim("PM2sources.txt", stringsAsFactors = F, sep="\t") %>%
  select(array, position, carbon_source)
PM3sources <- read.delim("PM3sources.txt", stringsAsFactors = F, sep="\t") %>%
  select(array, position, carbon_source)

# list the data files that were measured at 750nm.  The files are names as strain, time point, plate (cbs573t0pm1_750.txt)
PM1_r1_file_list <- list.files(path="../01_clean_files2/750/Dec2022/", pattern = glob2rx("*pm1*_750.txt"))
PM1_r2_file_list <- list.files(path="../01_clean_files2/750/June2023/", pattern = glob2rx("*pm1*_750.txt"))
PM2_r1_file_list <- list.files(path="../01_clean_files2/750/Dec2022/", pattern = glob2rx("*pm2*_750.txt"))
PM2_r2_file_list <- list.files(path="../01_clean_files2/750/June2023/", pattern = glob2rx("*pm2*_750.txt"))
PM3_r1_file_list <- list.files(path="../01_clean_files2/750/Dec2022/", pattern = glob2rx("*pm3*_750.txt"))
PM3_r2_file_list <- list.files(path="../01_clean_files2/750/June2023/", pattern = glob2rx("*pm3*_750.txt"))

PM1_euc_file_list_a <- list.files(path="../01_clean_files3/750/", pattern = glob2rx("*pm1*_750.txt"))
PM1_euc_file_list_b <- list.files(path="../01_clean_files4/620/", pattern = glob2rx("*pm1*_620.txt"))
PM2_euc_file_list_a <- list.files(path="../01_clean_files3/750/", pattern = glob2rx("*pm2*_750.txt"))
PM2_euc_file_list_b <- list.files(path="../01_clean_files4/620/", pattern = glob2rx("*pm2*_620.txt"))


# use a loop to combine all PM1 data
PM1_r1_data <- data.frame()
for (i in 1:length(PM1_r1_file_list)){
  temp_data <- read.delim(paste("~/Documents/2022_10_phenotypes/01_clean_files2/750/Dec2022/", PM1_r1_file_list[i], sep = ""), stringsAsFactors = F, header = FALSE)
  colnames(temp_data) = c("abs")
  temp_data <- temp_data %>%
    add_column(plate = PM1_r1_file_list[i], .before = "abs") %>%
    dplyr::mutate(plate = str_replace(plate, "(\\S*)(.txt)", "\\1")) %>%
    dplyr::relocate(plate)
  temp_data <- bind_cols(PM1sources, temp_data)
  PM1_r1_data <- rbindlist(list(PM1_r1_data, temp_data), use.names = T)
}

PM1_r2_data <- data.frame()
for (i in 1:length(PM1_r2_file_list)){
  temp_data <- read.delim(paste("~/Documents/2022_10_phenotypes/01_clean_files2/750/June2023/", PM1_r2_file_list[i], sep = ""), stringsAsFactors = F, header = FALSE)
  colnames(temp_data) = c("abs")
  temp_data <- temp_data %>%
    add_column(plate = PM1_r2_file_list[i], .before = "abs") %>%
    dplyr::mutate(plate = str_replace(plate, "(\\S*)(.txt)", "\\1")) %>%
    dplyr::relocate(plate)
  temp_data <- bind_cols(PM1sources, temp_data)
  PM1_r2_data <- rbindlist(list(PM1_r2_data, temp_data), use.names = T)
}

PM1_euc_data_a <- data.frame()
for (i in 1:length(PM1_euc_file_list_a)){
  temp_data <- read.delim(paste("~/Documents/2022_10_phenotypes/01_clean_files3/750/", PM1_euc_file_list_a[i], sep = ""), stringsAsFactors = F, header = FALSE)
  colnames(temp_data) = c("abs")
  temp_data <- temp_data %>%
    add_column(plate = PM1_euc_file_list_a[i], .before = "abs") %>%
    dplyr::mutate(plate = str_replace(plate, "(\\S*)(.txt)", "\\1")) %>%
    dplyr::relocate(plate)
  temp_data <- bind_cols(PM1sources, temp_data)
  PM1_euc_data_a <- rbindlist(list(PM1_euc_data_a, temp_data), use.names = T)
}

PM2_euc_data_a <- data.frame()
for (i in 1:length(PM2_euc_file_list_a)){
  temp_data <- read.delim(paste("~/Documents/2022_10_phenotypes/01_clean_files3/750/", PM2_euc_file_list_a[i], sep = ""), stringsAsFactors = F, header = FALSE)
  colnames(temp_data) = c("abs")
  temp_data <- temp_data %>%
    add_column(plate = PM2_euc_file_list_a[i], .before = "abs") %>%
    dplyr::mutate(plate = str_replace(plate, "(\\S*)(.txt)", "\\1")) %>%
    dplyr::relocate(plate)
  temp_data <- bind_cols(PM2sources, temp_data)
  PM2_euc_data_a <- rbindlist(list(PM2_euc_data_a, temp_data), use.names = T)
}

PM1_euc_data_b <- data.frame()
for (i in 1:length(PM1_euc_file_list_b)){
  temp_data <- read.delim(paste("~/Documents/2022_10_phenotypes/01_clean_files4/620/", PM1_euc_file_list_b[i], sep = ""), stringsAsFactors = F, header = FALSE)
  colnames(temp_data) = c("abs")
  temp_data <- temp_data %>%
    add_column(plate = PM1_euc_file_list_b[i], .before = "abs") %>%
    dplyr::mutate(plate = str_replace(plate, "(\\S*)(.txt)", "\\1")) %>%
    dplyr::relocate(plate)
  temp_data <- bind_cols(PM1sources, temp_data)
  PM1_euc_data_b <- rbindlist(list(PM1_euc_data_b, temp_data), use.names = T)
}

PM2_euc_data_b <- data.frame()
for (i in 1:length(PM2_euc_file_list_b)){
  temp_data <- read.delim(paste("~/Documents/2022_10_phenotypes/01_clean_files4/620/", PM2_euc_file_list_b[i], sep = ""), stringsAsFactors = F, header = FALSE)
  colnames(temp_data) = c("abs")
  temp_data <- temp_data %>%
    add_column(plate = PM2_euc_file_list_b[i], .before = "abs") %>%
    dplyr::mutate(plate = str_replace(plate, "(\\S*)(.txt)", "\\1")) %>%
    dplyr::relocate(plate)
  temp_data <- bind_cols(PM2sources, temp_data)
  PM2_euc_data_b <- rbindlist(list(PM2_euc_data_b, temp_data), use.names = T)
}

# use a loop to combine all PM2 data
PM2_r1_data <- data.frame()
for (i in 1:length(PM2_r1_file_list)){
  temp_data <- read.delim(paste("~/Documents/2022_10_phenotypes/01_clean_files2/750/Dec2022/", PM2_r1_file_list[i], sep = ""), stringsAsFactors = F, header = FALSE)
  colnames(temp_data) = c("abs")
  temp_data <- temp_data %>%
    add_column(plate = PM2_r1_file_list[i], .before = "abs") %>%
    dplyr::mutate(plate = str_replace(plate, "(\\S*)(.txt)", "\\1")) %>%
    dplyr::relocate(plate)
  temp_data <- bind_cols(PM2sources, temp_data)
  PM2_r1_data <- rbindlist(list(PM2_r1_data, temp_data), use.names = T)
}

PM2_r2_data <- data.frame()
for (i in 1:length(PM2_r2_file_list)){
  temp_data <- read.delim(paste("~/Documents/2022_10_phenotypes/01_clean_files2/750/June2023/", PM2_r2_file_list[i], sep = ""), stringsAsFactors = F, header = FALSE)
  colnames(temp_data) = c("abs")
  temp_data <- temp_data %>%
    add_column(plate = PM2_r2_file_list[i], .before = "abs") %>%
    dplyr::mutate(plate = str_replace(plate, "(\\S*)(.txt)", "\\1")) %>%
    dplyr::relocate(plate)
  temp_data <- bind_cols(PM2sources, temp_data)
  PM2_r2_data <- rbindlist(list(PM2_r2_data, temp_data), use.names = T)
}

# use a loop to combine all PM3 data
PM3_r1_data <- data.frame()
for (i in 1:length(PM3_r1_file_list)){
  temp_data <- read.delim(paste("~/Documents/2022_10_phenotypes/01_clean_files2/750/Dec2022/", PM3_r1_file_list[i], sep = ""), stringsAsFactors = F, header = FALSE)
  colnames(temp_data) = c("abs")
  temp_data <- temp_data %>%
    add_column(plate = PM3_r1_file_list[i], .before = "abs") %>%
    dplyr::mutate(plate = str_replace(plate, "(\\S*)(.txt)", "\\1")) %>%
    dplyr::relocate(plate)
  temp_data <- bind_cols(PM3sources, temp_data)
  PM3_r1_data <- rbindlist(list(PM3_r1_data, temp_data), use.names = T)
}

PM3_r2_data <- data.frame()
for (i in 1:length(PM3_r2_file_list)){
  temp_data <- read.delim(paste("~/Documents/2022_10_phenotypes/01_clean_files2/750/June2023/", PM3_r2_file_list[i], sep = ""), stringsAsFactors = F, header = FALSE)
  colnames(temp_data) = c("abs")
  temp_data <- temp_data %>%
    add_column(plate = PM3_r2_file_list[i], .before = "abs") %>%
    dplyr::mutate(plate = str_replace(plate, "(\\S*)(.txt)", "\\1")) %>%
    dplyr::relocate(plate)
  temp_data <- bind_cols(PM3sources, temp_data)
  PM3_r2_data <- rbindlist(list(PM3_r2_data, temp_data), use.names = T)
}

# combine the arrays
PM_r1_data <- rbindlist(list(PM1_r1_data, PM2_r1_data, PM3_r1_data), use.names = T)
PM_r2_data <- rbindlist(list(PM1_r2_data, PM2_r2_data, PM3_r2_data), use.names = T)
PM_euc_data <- rbindlist(list(PM1_euc_data_a, PM2_euc_data_a), use.names = T)

# put strain information in its own column.  Add a timepoint column "T".
PM_r1_data <- PM_r1_data %>%
  dplyr::mutate(organism = str_replace(plate, "(\\S*)(t\\d)(pm\\d)(r\\d)(_750)", "\\1")) %>%
  dplyr::mutate(T = str_replace(plate, "(\\S*)(t)(\\d)(pm\\d)(r\\d)(_750)", "\\3")) %>%
  dplyr::mutate(T = as.numeric(T)) %>%
  dplyr::mutate(Rep = str_replace(plate, "(\\S*)(t)(\\d)(pm\\d)(r)(\\d)(_750)", "\\6")) %>%
  dplyr::mutate(Rep = as.numeric(Rep))
str(PM_r1_data)

PM_r2_data <- PM_r2_data %>%
  dplyr::mutate(organism = str_replace(plate, "(\\S*)(t\\d)(pm\\d)(r\\d)(_750)", "\\1")) %>%
  dplyr::mutate(T = str_replace(plate, "(\\S*)(t)(\\d)(pm\\d)(r\\d)(_750)", "\\3")) %>%
  dplyr::mutate(T = as.numeric(T)) %>%
  dplyr::mutate(Rep = str_replace(plate, "(\\S*)(t)(\\d)(pm\\d)(r)(\\d)(_750)", "\\6")) %>%
  dplyr::mutate(Rep = as.numeric(Rep))

PM_euc_data <- PM_euc_data %>%
  dplyr::mutate(organism = str_replace(plate, "(\\S*)(t\\d)(pm\\d)(r\\d)(_\\d\\d\\d)", "\\1")) %>%
  dplyr::mutate(T = str_replace(plate, "(\\S*)(t)(\\d)(pm\\d)(r\\d)(_\\d\\d\\d)", "\\3")) %>%
  dplyr::mutate(T = as.numeric(T)) %>%
  dplyr::mutate(Rep = str_replace(plate, "(\\S*)(t)(\\d)(pm\\d)(r)(\\d)(_\\d\\d\\d)", "\\6")) %>%
  dplyr::mutate(Rep = as.numeric(Rep))

# add in the time in hours (as opposed to numbered timepoints)
time_points_r1 = read.delim("time_points_r1.txt", stringsAsFactors = F, sep="\t")
PM_r1_data = dplyr::left_join(PM_r1_data, time_points_r1, by = c("T"))

time_points_r2 = read.delim("time_points_r2.txt", stringsAsFactors = F, sep="\t")
PM_r2_data = dplyr::left_join(PM_r2_data, time_points_r2, by = c("T"))

time_points_euc = read.delim("time_points_euc.txt", stringsAsFactors = F, sep="\t")
PM_euc_data = dplyr::left_join(PM_euc_data, time_points_euc, by = c("T"))

PM_data <- rbindlist(list(PM_r1_data, PM_r2_data, PM_euc_data), use.names = T)

str(PM_data)

carb_data <- PM_data %>% filter(array != "PM3B")
length(unique(carb_data$carbon_source))

# checking that I have the right number of elements in PM3 (96)
nitro_data <- PM_data %>% filter(array == "PM3B")
length(unique(nitro_data$carbon_source))

########################


# make a named vector of viridis colours named by carbon source
carbon_source <- unique(carb_data$carbon_source) %>% sort()
carbon_source <- carbon_source[ !carbon_source == "None"]
carbon_source
colours_carb <- viridis::viridis(n = length(carbon_source))
names(colours_carb) <- carbon_source
colours_carb
str(colours_carb)
negative <- c("None"="black")
negative
str(negative)

colours_carb <- c(colours_carb, negative)


# Make a function to plot it out ------------------------------------------


# this replaces the first section, so you can just run this instead of all that code
plot_carbon <- function(df, target_array, target_organism, target_rep, target_cutoff, colour = colours_carb){

  df_1 <- df %>% 
    filter(array == target_array, organism == target_organism, Rep == target_rep) %>% 
    group_by(carbon_source) %>% 
    mutate(over_20 = if_else(any(hours == max(hours) & abs > paste(target_cutoff)), 'yes', 'no')) %>% 
    ungroup()
  df_1 <- df_1 %>% mutate(over_20 = ifelse(carbon_source == "None" & over_20 == "no", "yes", over_20))
  
  annotations <- df_1 %>% 
    filter(hours == max(hours), over_20 == 'yes') %>% 
    mutate(hours = 169)
  
  #negative <- df_1 %>% 
   # filter(hours == max(hours), carbon_source == 'Negative Control') %>% 
  #  mutate(hours = 169)
  
  #annotations <- rbind(annotations_1, negative)
  
  p <- df_1 %>% 
    ggplot(aes(x = hours, y = abs, col = carbon_source)) +
    #geom_hline(yintercept = paste(target_cutoff), linetype = 'dashed', alpha = 0.5) +
    ggtitle(paste(target_array," ",target_organism, " ", "rep", target_rep, sep = "")) +
    geom_line(aes(alpha = over_20)) +
    ggrepel::geom_text_repel(data = annotations, 
                             aes(label = carbon_source),
                             hjust = -0.1,
                             direction = 'y',
                             min.segment.length = Inf) + 
    xlim(c(0, 300)) + 
    scale_alpha_discrete(range = c(0.2, 1)) + 
    scale_color_manual(values = colours_carb) + 
    theme(plot.title = element_text(color="black", size=18, face="bold", hjust = 0.5),
          legend.position = "none",
          axis.text.y = element_text(colour = "black", size = 12), 
          axis.text.x = element_text(colour = "black", size = 12),
          axis.title.y = element_text(size = 18), 
          axis.title.x = element_text(size = 18, colour = "black"),
          panel.background = element_blank())
  
  return(p)
}

# run plots

PM1_CBS119000_r1 = plot_carbon(df = PM_data, target_array = 'PM1', target_organism = 'cbs119000', target_rep = '1', target_cutoff = '0.375')
PM1_CBS119000_r1
PM1_CBS119000_r2 = plot_carbon(df = PM_data, target_array = 'PM1', target_organism = 'cbs119000', target_rep = '2', target_cutoff = '0.3')
PM1_CBS119000_r2
PM1_CBS119000_r3 = plot_carbon(df = PM_data, target_array = 'PM1', target_organism = 'cbs119000', target_rep = '3', target_cutoff = '0.2')
PM1_CBS119000_r3

PM1_CBS57363_r1 = plot_carbon(df = PM_data, target_array = 'PM1', target_organism = 'cbs57363', target_rep = '1', target_cutoff = '0.375')
PM1_CBS57363_r1
PM1_CBS57363_r2 = plot_carbon(df = PM_data, target_array = 'PM1', target_organism = 'cbs57363', target_rep = '2', target_cutoff = '0.3')
PM1_CBS57363_r2
PM1_CBS57363_r3 = plot_carbon(df = PM_data, target_array = 'PM1', target_organism = 'cbs57363', target_rep = '3', target_cutoff = '0.2')
PM1_CBS57363_r3

PM1_semc78_r1 = plot_carbon(df = PM_data, target_array = 'PM1', target_organism = 'semc78', target_rep = '1', target_cutoff = '0.2')
PM1_semc78_r1
PM1_semc78_r2 = plot_carbon(df = PM_data, target_array = 'PM1', target_organism = 'semc78', target_rep = '2', target_cutoff = '0.2')
PM1_semc78_r2
PM1_semc78_r3 = plot_carbon(df = PM_data, target_array = 'PM1', target_organism = 'semc78', target_rep = '3', target_cutoff = '0.2')
PM1_semc78_r3

PM1_semc49_r1 = plot_carbon(df = PM_data, target_array = 'PM1', target_organism = 'semc49', target_rep = '1', target_cutoff = '0.2')
PM1_semc49_r1
PM1_semc49_r2 = plot_carbon(df = PM_data, target_array = 'PM1', target_organism = 'semc49', target_rep = '2', target_cutoff = '0.2')
PM1_semc49_r2
PM1_semc49_r3 = plot_carbon(df = PM_data, target_array = 'PM1', target_organism = 'semc49', target_rep = '3', target_cutoff = '0.2')
PM1_semc49_r3

PM1_euc_r1 = plot_carbon(df = PM_data, target_array = 'PM1', target_organism = 'cbs140593', target_rep = '1', target_cutoff = '0.4')
PM1_euc_r1
PM1_euc_r2 = plot_carbon(df = PM_data, target_array = 'PM1', target_organism = 'cbs140593', target_rep = '2', target_cutoff = '0.4')
PM1_euc_r2
PM1_euc_r3 = plot_carbon(df = PM_data, target_array = 'PM1', target_organism = 'cbs140593', target_rep = '3', target_cutoff = '0.4')
PM1_euc_r3

PM2_CBS119000_r1 = plot_carbon(df = PM_data, target_array = 'PM2A', target_organism = 'cbs119000', target_rep = '1', target_cutoff = '0.2')
PM2_CBS119000_r1
PM2_CBS119000_r2 = plot_carbon(df = PM_data, target_array = 'PM2A', target_organism = 'cbs119000', target_rep = '2', target_cutoff = '0.2')
PM2_CBS119000_r2
PM2_CBS119000_r3 = plot_carbon(df = PM_data, target_array = 'PM2A', target_organism = 'cbs119000', target_rep = '3', target_cutoff = '0.2')
PM2_CBS119000_r3

PM2_CBS57363_r1 = plot_carbon(df = PM_data, target_array = 'PM2A', target_organism = 'cbs57363', target_rep = '1', target_cutoff = '0.2')
PM2_CBS57363_r1
PM2_CBS57363_r2 = plot_carbon(df = PM_data, target_array = 'PM2A', target_organism = 'cbs57363', target_rep = '2', target_cutoff = '0.2')
PM2_CBS57363_r2
PM2_CBS57363_r3 = plot_carbon(df = PM_data, target_array = 'PM2A', target_organism = 'cbs57363', target_rep = '3', target_cutoff = '0.2')
PM2_CBS57363_r3

PM2_semc78_r1 = plot_carbon(df = PM_data, target_array = 'PM2A', target_organism = 'semc78', target_rep = '1', target_cutoff = '0.3')
PM2_semc78_r1
PM2_semc78_r2 = plot_carbon(df = PM_data, target_array = 'PM2A', target_organism = 'semc78', target_rep = '2', target_cutoff = '0.3')
PM2_semc78_r2
PM2_semc78_r3 = plot_carbon(df = PM_data, target_array = 'PM2A', target_organism = 'semc78', target_rep = '3', target_cutoff = '0.3')
PM2_semc78_r3

PM2_semc49_r1 = plot_carbon(df = PM_data, target_array = 'PM2A', target_organism = 'semc49', target_rep = '1', target_cutoff = '0.2')
PM2_semc49_r1
PM2_semc49_r2 = plot_carbon(df = PM_data, target_array = 'PM2A', target_organism = 'semc49', target_rep = '2', target_cutoff = '0.2')
PM2_semc49_r2
PM2_semc49_r3 = plot_carbon(df = PM_data, target_array = 'PM2A', target_organism = 'semc49', target_rep = '3', target_cutoff = '0.2')
PM2_semc49_r3

PM2_euc_r1 = plot_carbon(df = PM_data, target_array = 'PM2A', target_organism = 'cbs140593', target_rep = '1', target_cutoff = '0.3')
PM2_euc_r1
PM2_euc_r2 = plot_carbon(df = PM_data, target_array = 'PM2A', target_organism = 'cbs140593', target_rep = '2', target_cutoff = '0.3')
PM2_euc_r2
PM2_euc_r3 = plot_carbon(df = PM_data, target_array = 'PM2A', target_organism = 'cbs140593', target_rep = '3', target_cutoff = '0.3')
PM2_euc_r3



# make a composite plot

library(cowplot)
pdf(file="three_reps/PM1_epigloea.pdf", width=20, height=10)
plot_grid(PM1_CBS119000_r1, PM1_CBS119000_r2, PM1_CBS119000_r3, PM1_CBS57363_r1, PM1_CBS57363_r2, PM1_CBS57363_r3, ncol = 3)
dev.off()

pdf(file="three_reps/PM1_eucalyptigena.pdf", width=20, height=10)
plot_grid(PM1_euc_r1, PM1_euc_r2, PM1_euc_r3, ncol = 3)
dev.off()

pdf(file="three_reps/PM1_Sporothrix.pdf", width=20, height=10)
plot_grid(PM1_CBS119000_r1, PM1_CBS119000_r2, PM1_CBS119000_r3, PM1_CBS57363_r1, PM1_CBS57363_r2, PM1_CBS57363_r3, PM1_euc_r1, PM1_euc_r2, PM1_euc_r3, ncol = 3)
dev.off()

pdf(file="three_reps/PM1_associates.pdf", width=20, height=10)
plot_grid(PM1_semc49_r1, PM1_semc49_r2, PM1_semc49_r3, PM1_semc78_r1, PM1_semc78_r2, PM1_semc78_r3, ncol = 3)
dev.off()

pdf(file="three_reps/PM1_all.pdf", width=20, height=20)
plot_grid(PM1_CBS119000_r1, PM1_CBS119000_r2, PM1_CBS119000_r3, PM1_CBS57363_r1, PM1_CBS57363_r2, PM1_CBS57363_r3, PM1_semc49_r1, PM1_semc49_r2, PM1_semc49_r3, PM1_semc78_r1, PM1_semc78_r2, PM1_semc78_r3, ncol = 3)
dev.off()

pdf(file="three_reps/PM2_epigloea.pdf", width=20, height=10)
plot_grid(PM2_CBS119000_r1, PM2_CBS119000_r2, PM2_CBS119000_r3, PM2_CBS57363_r1, PM2_CBS57363_r2, PM2_CBS57363_r3, ncol = 3)
dev.off()

pdf(file="three_reps/PM2_eucalyptigena.pdf", width=20, height=10)
plot_grid(PM2_euc_r1, PM2_euc_r2, PM2_euc_r3, ncol = 3)
dev.off()

pdf(file="three_reps/PM2_Sporothrix.pdf", width=20, height=10)
plot_grid(PM2_CBS119000_r1, PM2_CBS119000_r2, PM2_CBS119000_r3, PM2_CBS57363_r1, PM2_CBS57363_r2, PM2_CBS57363_r3, PM2_euc_r1, PM2_euc_r2, PM2_euc_r3, ncol = 3)
dev.off()

pdf(file="three_reps/PM2_associates.pdf", width=20, height=10)
plot_grid(PM2_semc78_r1, PM2_semc78_r2, PM2_semc78_r3, PM2_semc49_r1, PM2_semc49_r2, PM2_semc49_r3, ncol = 3)
dev.off()

pdf(file="three_reps/PM2_all.pdf", width=20, height=20)
plot_grid(PM2_CBS119000_r1, PM2_CBS119000_r2, PM2_CBS119000_r3, PM2_CBS57363_r1, PM2_CBS57363_r2, PM2_CBS57363_r3, PM2_semc49_r1, PM2_semc49_r2, PM2_semc49_r3, PM2_semc78_r1, PM2_semc78_r2, PM2_semc78_r3, ncol = 3)
dev.off()



############################################

# Nitrogen array

# make a named vector of viridis colours named by nitrogen source
nitrogen_source <- unique(nitro_data$carbon_source) %>% sort()
nitrogen_source <- nitrogen_source[ !nitrogen_source == "Negative Control"]
colours_nitro <- viridis::viridis(n = length(nitrogen_source))
names(colours_nitro) <- nitrogen_source
colours_nitro
str(colours_nitro)
negative <- c("Negative Control"="black")
negative
str(negative)

colours_nitro <- c(colours_nitro, negative)


# this replaces the first section, so you can just run this instead of all that code
plot_nitrogen <- function(df, target_array, target_organism, target_rep, target_cutoff, colour = colours_nitro){
  
  df_1 <- df %>% 
    filter(array == target_array, organism == target_organism, Rep == target_rep) %>% 
    group_by(carbon_source) %>% 
    mutate(over_20 = if_else(any(hours == max(hours) & abs > paste(target_cutoff)), 'yes', 'no')) %>% 
    ungroup()
  df_1 <- df_1 %>% mutate(over_20 = ifelse(carbon_source == "Negative Control" & over_20 == "no", "yes", over_20))
  
  annotations <- df_1 %>% 
    filter(hours == max(hours), over_20 == 'yes') %>% 
    mutate(hours = 169)
  
  #negative <- df_1 %>% 
  # filter(hours == max(hours), carbon_source == 'Negative Control') %>% 
  #  mutate(hours = 169)
  
  #annotations <- rbind(annotations_1, negative)
  
  p <- df_1 %>% 
    ggplot(aes(x = hours, y = abs, col = carbon_source)) +
    #geom_hline(yintercept = paste(target_cutoff), linetype = 'dashed', alpha = 0.5) +
    ggtitle(paste(target_array," ",target_organism, " ", "rep", target_rep, sep = "")) +
    geom_line(aes(alpha = over_20)) +
    ggrepel::geom_text_repel(data = annotations, 
                             aes(label = carbon_source),
                             hjust = -0.1,
                             direction = 'y')+ 
    xlim(c(0, 300)) + 
    scale_alpha_discrete(range = c(0.2, 1)) + 
    scale_color_manual(values = colours_nitro) + 
    theme(plot.title = element_text(color="black", size=18, face="bold", hjust = 0.5),
          legend.position = "none",
          axis.text.y = element_text(colour = "black", size = 12), 
          axis.text.x = element_text(colour = "black", size = 12),
          axis.title.y = element_text(size = 18), 
          axis.title.x = element_text(size = 18, colour = "black"),
          panel.background = element_blank())
  
  return(p)
}

# run plots

PM3_CBS119000_r1 = plot_nitrogen(df = PM_data, target_array = 'PM3B', target_organism = 'cbs119000', target_rep = '1', target_cutoff = '0.3')
PM3_CBS119000_r1
PM3_CBS119000_r2 = plot_nitrogen(df = PM_data, target_array = 'PM3B', target_organism = 'cbs119000', target_rep = '2', target_cutoff = '0.3')
PM3_CBS119000_r2
PM3_CBS119000_r3 = plot_nitrogen(df = PM_data, target_array = 'PM3B', target_organism = 'cbs119000', target_rep = '3', target_cutoff = '0.3')
PM3_CBS119000_r3

PM3_CBS57363_r1 = plot_nitrogen(df = PM_data, target_array = 'PM3B', target_organism = 'cbs57363', target_rep = '1', target_cutoff = '0.3')
PM3_CBS57363_r1
PM3_CBS57363_r2 = plot_nitrogen(df = PM_data, target_array = 'PM3B', target_organism = 'cbs57363', target_rep = '2', target_cutoff = '0.3')
PM3_CBS57363_r2
PM3_CBS57363_r3 = plot_nitrogen(df = PM_data, target_array = 'PM3B', target_organism = 'cbs57363', target_rep = '3', target_cutoff = '0.3')
PM3_CBS57363_r3

PM3_semc78_r1 = plot_nitrogen(df = PM_data, target_array = 'PM3B', target_organism = 'semc78', target_rep = '1', target_cutoff = '0.3')
PM3_semc78_r1
PM3_semc78_r2 = plot_nitrogen(df = PM_data, target_array = 'PM3B', target_organism = 'semc78', target_rep = '2', target_cutoff = '0.3')
PM3_semc78_r2
PM3_semc78_r3 = plot_nitrogen(df = PM_data, target_array = 'PM3B', target_organism = 'semc78', target_rep = '3', target_cutoff = '0.3')
PM3_semc78_r3

PM3_semc49_r1 = plot_nitrogen(df = PM_data, target_array = 'PM3B', target_organism = 'semc49', target_rep = '1', target_cutoff = '0.3')
PM3_semc49_r1
PM3_semc49_r2 = plot_nitrogen(df = PM_data, target_array = 'PM3B', target_organism = 'semc49', target_rep = '2', target_cutoff = '0.3')
PM3_semc49_r2
PM3_semc49_r3 = plot_nitrogen(df = PM_data, target_array = 'PM3B', target_organism = 'semc49', target_rep = '3', target_cutoff = '0.3')
PM3_semc49_r3

# make a composite plot

library(cowplot)
pdf(file="three_reps/PM3_epigloea.pdf", width=20, height=10)
plot_grid(PM3_CBS119000_r1, PM3_CBS119000_r2, PM3_CBS119000_r3, PM3_CBS57363_r1, PM3_CBS57363_r2, PM3_CBS57363_r3, ncol = 3)
dev.off()

pdf(file="three_reps/PM3_associates.pdf", width=20, height=10)
plot_grid(PM3_semc49_r1, PM3_semc49_r2, PM3_semc49_r3, PM3_semc78_r1, PM3_semc78_r2, PM3_semc78_r3, ncol = 3)
dev.off()


pdf(file="three_reps/PM3_all.pdf", width=40, height=40)
plot_grid(PM3_CBS119000_r1, PM3_CBS119000_r2, PM3_CBS119000_r3, PM3_CBS57363_r1, PM3_CBS57363_r2, PM3_CBS57363_r3, PM3_semc49_r1, PM3_semc49_r2, PM3_semc49_r3, PM3_semc78_r1, PM3_semc78_r2, PM3_semc78_r3, ncol = 3)
dev.off()


sessionInfo()
