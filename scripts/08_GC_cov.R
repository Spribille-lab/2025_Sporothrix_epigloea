#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

sample <- args[1]
sample <- "CBS95973"

#This script produces t GC/coverage plots for a single metagenome.
#As input, it uses for each metagenome: 
#1.nucleotide assembly (contigs.fasta`)
#2. CONCOCT-produced coverage files (coverage_table.tsv)
#3. CONCOCT-produced binning files (clustering_gt1000_merged.csv)

setwd(paste("~/Documents/2023_02_Sporothrix/08_GC_coverage/",sample,"/", sep=""))

getwd()


library(tidyverse)
library(Biostrings)
library(DT)

#create functions
getGC<-function(fasta){
  #process fasta
  contig_id<-names(fasta) #get contig names
  nucl_freq<-alphabetFrequency(fasta)
  nucl_freq<-as.data.frame(nucl_freq) %>% mutate(gc=100*(G+C+S)/(A+T+W+G+C+S)) #get gc content
  dataset1<-data.frame(contig_id,nucl_freq$gc)
  return(dataset1)
}

####DATA PROCESSING
#load files
depth1<-read.delim("data/coverage_table.tsv")
cluster<-read.csv("data/clustering_gt1000_merged.csv")
fasta<-readDNAStringSet('data/contigs.fasta')

#process depth file
colnames(depth1)<-c("contig_region","coverage")
depth1<-transform(depth1, contig_id = sub("(.*_cov_[0-9]+\\.[0-9]+).*", "\\1", contig_region)) #depth file gives coverage for small chunks of a contig. the name of chunks follow the scheme "contig_id.1", "contig_id.2", etc. This line uses regular expression to extract the contig_id from each chunk name
depth2<-depth1 %>% group_by(contig_id) %>% summarize(mean_contig_cov=mean(coverage)) #calculates mean coverage for each contig

#combine depth and binning
df<-left_join(depth2,cluster)

#add GC content
gc<-getGC(fasta)
df2<-left_join(df,gc)

#add contig length
length<-data.frame(names(fasta),width(fasta))
colnames(length)<-c("contig_id","length")
df2<-left_join(df2,length)

#save the final table
write.table(df2,"reports/gc_cov.txt",sep='\t',quote = F,row.names = F)

#calculate and save mean coverage per bin
bin_coverage<-df2 %>% group_by(cluster_id)%>%summarize(cov=mean(mean_contig_cov))
write.table(bin_coverage,"reports/mean_bin_coverage.txt",sep='\t',quote = F,row.names = F)

#calculate and save mean GC per bin
bin_gc<-df2 %>% group_by(cluster_id)%>%summarize(gc=mean(nucl_freq.gc))

# combine mean coverage and GC
df3<-left_join(bin_coverage,bin_gc)

# make gc_cov plot with contigs colored according to their bins

library(plotly)
gc_plot<-plot_ly(data=df2,x=~nucl_freq.gc,y=~mean_contig_cov,color=~as.factor(cluster_id),
                 marker=list(size=4)) %>%
  add_markers() %>%  
  layout( xaxis = list(title = 'GC%'), 
          yaxis = list(title = 'Coverage',type="log"))

#save as html
withr::with_dir('reports', htmlwidgets::saveWidget(as_widget(gc_plot), file="gc_plot.html"))
#save as png - optional requires installing orca (https://github.com/plotly/orca)
#Sys.setenv("PATH" = paste(Sys.getenv("PATH"), "/opt/anaconda2/bin/", sep = .Platform$path.sep))
#orca(gc_plot,"figures/gc_plot.svg")


