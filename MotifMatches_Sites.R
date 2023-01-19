#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Enter path to MotifMatches table", call.=FALSE)
} 

###Libraries

library(tidyverse)
library(gridExtra)

#Read MotifMatches_Dis.py results table
#input <- 'MotifMatches_list.txt'
input1 <- args[1]
results_table <- gsub('_list.txt','',input1)
dat <- read.table(input1, header=T,sep="\t")

#Create a list of the motif matches sites to retrieve alignment presence 
dat_sites <- dat[,c(1:3,5)]
dat_sites$Motif_sSite<-as.character(dat_sites$Motif_sSite)
dat_sites<- dat_sites %>% 
  group_by(Protein_ID,Motif_Name)  %>% 
  mutate(Motif_sites = paste(Motif_sSite, collapse = ","))
dat_sites<-select(dat_sites,-4)
dat_sites<-dat_sites[!duplicated(dat_sites),]

write_tsv(dat_sites, file=paste(results_table,'_sites.txt',sep = ""),col_names = T)

rm(dat,dat_sites)
rm(input1)
