#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Enter path to MotifMatches table & LOPIT file ", call.=FALSE)
} 

###Libraries

library(tidyverse)
library(gridExtra)

#Define working directory
setwd("/Users/JAVlvrd/Documents/Toxoplasma-2022/ToxoMotifs/Results/")

#Alias for reading further files
input1<- args[1]
# input2<- args[2]
# input3<- args[3]
# input4<- args[4]
#input1 <- "ELM_test"

#Motif Matches file processing for data integration
#input1 <- args[1]
#input1 <- "/Users/JAVlvrd/Documents/Toxoplasma-2022/ToxoMotifs/Results/ELM_Sep22_MotifMatches_list.txt"
matches_file <- paste(input1,"_MotifMatches_list.txt",sep="")
#matches_file <- input1
motif_list<- read_tsv(matches_file)
motif_list$seq_id<-gsub("-t26_1-p1","",motif_list$Protein_ID ) #Define common sequence ID key for further data integration
motif_list$key <- paste(motif_list$seq_id,motif_list$Motif_Name,motif_list$Match_N,sep="|") #key should be a combination of sequence ID, motif name and motif instance
motif_list$Motif_Type <- substr(motif_list$Motif_Name,1,3)

##### Add MSA scores ####
#Read conservation presence file and combine it with Motif list
cons_file <- paste(input1,"_MotifMatches_presence.txt",sep="")
#cons_file <- input2
cons_table <- read_tsv(cons_file,col_names =T)
cons_table <- select(cons_table, -c(2:4))
rm(cons_file)

MotifHits <-  left_join(motif_list,cons_table, by = "key") 
rm(cons_table)

#Add Expression evidence ##Alingments were only done with proteins with expression evidence MassSpecEv Peptides>=1
MotifHits$ExpEvd <- sapply(MotifHits$presence_org,function(x) if(is.na(x)){'no'}else{"yes"})

##### Add LOPIT locations ####
#Read  LOPIT location  table
Tgondii_LOPIT <- read_tsv(file = '/Users/JAVlvrd/Documents/Toxoplasma-2021/AminoAcidComposition/RefTables/TgME49_LOPIT_MS_evidence.txt') 
#TAGM-MAP (default probability) > 0.4 % 07.06.22 MassSpecEv Peptides>=1
Tgondii_LOPIT <- Tgondii_LOPIT %>% rename(seq_id = Gene_ID) 
Tgondii_LOPIT <- select(Tgondii_LOPIT, -c(2,6:7))
MotifHits <- left_join(MotifHits,Tgondii_LOPIT, by = "seq_id") 
rm(Tgondii_LOPIT)

#Add Location evidence
MotifHits$LocEvd <- sapply(MotifHits$Predicted_Location_TAGM_MAP,function(x) if(is.na(x)){'no'}else{"yes"})

#Add a summary of locations
MotifHits$Organelle <- gsub("_[1-9]","",MotifHits$Predicted_Location_TAGM_MAP)
MotifHits$Organelle <- gsub("_-_.*?$","",MotifHits$Organelle)
MotifHits$Organelle <- gsub("^[0-9][0-9]S_","",MotifHits$Organelle)


#### Add Bradyzoite info ######
#https://doi.org/10.1371/journal.pone.0232552
brad <- read_tsv("/Users/JAVlvrd/Documents/Toxoplasma-2021/Motif_Enrichments/Nadipuram2020_GRAnames.txt")
tach <- MotifHits %>% filter(Organelle=="dense_granules") %>% select(seq_id) %>% unique()
temp <- intersect(tach$seq_id, brad$ToxoDB)
temp <- setdiff( brad$ToxoDB,tach$seq_id)

bradMotifs <- MotifHits[MotifHits$seq_id %in% temp,]
bradMotifs <- bradMotifs %>% select(seq_id,LocEvd,Organelle,Product_Description)
bradMotifs <- bradMotifs[!duplicated(bradMotifs$seq_id),]
#bradMotifs %>% select(Organelle) %>% table()

#update values with NA
#if the location is diferent it remains the one from the LOPIT
if(nrow(bradMotifs)>0){
  for(i in 1:length(temp)){
    change <-MotifHits[MotifHits$seq_id==temp[i],] %>% select(LocEvd)
    if(change[1,1]=="no"){
      MotifHits[MotifHits$seq_id==temp[i],19] <- "yes"
    }
    change <-MotifHits[MotifHits$seq_id==temp[i],]  %>% select(Organelle)
    if(is.na(change[1,1])){
      MotifHits[MotifHits$seq_id==temp[i],20] <- "dense_granules"
    }
    change <-MotifHits[MotifHits$seq_id==temp[i],]  %>% select(Product_Description)
    if(is.na(change[1,1])){
      MotifHits[MotifHits$seq_id==temp[i],15] <- brad[brad$ToxoDB==temp[i],2]
    }
  }
  rm(i, temp,change)
}
rm(bradMotifs, brad, tach)

#### Add phospho site info ####
##Collected from ToxoDB from PMID: 30850550, PMID: 21980283, PMID: 22018241
mods_file <- paste(input1,"_MotifMatches_mods.txt",sep="")
#mods_file <- input3
mods_table <-read_tsv(mods_file,col_names =T)
rm(mods_file)

MotifHits <-  left_join(MotifHits,mods_table, by = "key") 
rm(mods_table)


#### Add AlphaFold values ####
af_file <- paste(input1,"_MotifMatches_AF_extension.txt",sep="")
#af_file <- input4
af_table <- read_tsv(file=af_file, col_names = F)
colnames(af_table)<- c('key','mean_pLDDT','mean_acc')
af_table$mean_pLDDT<-as.numeric(af_table$mean_pLDDT)
af_table$mean_acc<-as.numeric(af_table$mean_acc)
rm(af_file)

MotifHits <- left_join(MotifHits,af_table, by="key")
rm(af_table)
# #### Add ToxoDB descriptions ####
# descriptions <- read_tsv("/Users/JAVlvrd/Documents/Toxoplasma-2021/Motif_Enrichments/result_Files/elm_classes_11.2021_filt_ToxoDB_descriptions.txt",col_names = F)
# descriptions <- select(descriptions, X1,X3)
# colnames(descriptions)<-c("seq_id","description")
# MotifHits <- left_join(MotifHits, descriptions, by="seq_id")

#### Save data #####
write_tsv(MotifHits, file=paste(input1,"_MotifMatches_complete.txt",sep=""),col_names = T)

#write_tsv(dat[,c(7,1:3,5,8,10)], file=paste(results_table,'_dat_short.txt',sep = ""),col_names = T)
#write_tsv(dat, file=paste(results_table,'_dat_full.txt',sep = ""),col_names = T)






