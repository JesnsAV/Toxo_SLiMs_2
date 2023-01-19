#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 16:52:39 2020

@author: JAVlvrd

It will take as input a Toxo Proteome and a list of motifs

It will give as output a list of motif hits with information about their disorder probabilities

"""


'''
PATH OF PROTEOME
'''

import argparse

parser=argparse.ArgumentParser()
parser.add_argument("data_alias",
                    help="Name alias for result files") #receive ALIAS
parser.add_argument("input_proteome",
                    help="Path to a proteome file") #receive proteome path
parser.add_argument("input_ELMs",
                    help="Path to a motif table file with names and REGEX") #receive file of signalP results
parser.add_argument("input_IUP_thrs", 
                    help="IUPRED pLDDT threshold", nargs='?', type=float, const=0.4) #receive file of signalP results
args = parser.parse_args()

#sample command
#python ../Software/MotifMatches_Dis.py 'ALIAS' ../Data/ToxoDB-42_TgondiiME49_AnnotatedProteins.fasta ../Data/elm_classes_200922.tsv 0.4

#outputs:
    #ALIAS_MotifMatches_list.txt #motif matches in a list
 
'''
Make a list of all proteins in the proteome
IDs are the dictionary keys
Information store include gene name, seq length and sequence
'''
print("Reading proteome...")

Proteome = args.input_proteome 
proteome_file = open(Proteome, 'r') #open up file
info_proteins = {}


for line in proteome_file:
    if line.startswith('>'): #identify the lines with header with '>'
        header = line.strip().split('|') #Separate the different pieces of information and store them in a list

        ID = header[0] #select the ID of the Sequence
        ID = ID.replace('>','').replace(' ','') #take away extra characters and spaces

        gene = header[2] #select the gene name of the Sequence
        gene = gene.replace(' gene=','').replace(' ','') #take away extra characters and spaces

        info_proteins[ID]=[gene, 0,''] #crate an entry in the info_proteins dic. using ID as key and storing gene name, aa length and empty sequence
        #adding an empty string to store amino acid sequence by extension
        del gene#, length #delete used variables not ID/key to keep track of which protein we are adding the seq_piece later

    else: #if the line is not a header then add it to the previous entry seq string
        seq_piece = line.strip() #take away extra formating from the line with aa's
        info_proteins[ID][2] = info_proteins[ID][2] + seq_piece #add currrent line to extend the protein sequence
        del seq_piece #delete
    info_proteins[ID][1]=len(info_proteins[ID][2])
del line, ID #delete used variables

keys = list(info_proteins.keys())

'''
Upload ELM models

'''
print("Loading ELM models...")

ELMs = args.input_ELMs 
ELM_table = open(ELMs, 'r') #open up file
ELM_models = {}


for line in ELM_table:
    row = line.strip().split() #Separate the name and the regular expression
    Name = row[0] #select the Name of the ELM
    REGEX = row[1] #select the REGEX of the ELM
    Group = row[2] #select the slim group

    ELM_models[Name]=[REGEX, Group]
    del row
del line, Name, Group #delete used variables

ELM_keys = list(ELM_models.keys()) #Get ELM dict keys



'''
Calculate Disorder of each aa position in the proteins provided
'''
print("Calculating protein disorder...")

import os
import numpy as np
import sys

sys.path.insert(1, '/Toxo_Proteomes') #path to python software
from iupred2a_folder.iupred2a import *

proteins_dis={} #create list that will contain all vectores with disorder values


ProteinList={}

for ID in keys:
    aa_seq = info_proteins[ID][2]  #save the aa sequence in a string variable
    aa_dis = iupred(aa_seq)[0] #calculate the disorder score (LONG default) for each aa in the sequence and save them in a list
    prot_dis = np.mean(aa_dis) #calculate the mean disorder of each protein and save it in a interger variable

    proteins_dis[ID]=aa_dis #save disorder score list in a dicctionary for later access
    ProteinList[ID]=prot_dis  #save average disorder in the general protein seq dictionary

    del aa_seq, aa_dis, prot_dis #delete variables after storing them



'''
Find motif matches and their average disorder
'''
print("Finding motif matches...")

import re

proteins_motifs={} #create an empty dict to save motif matches


match_id=0 #ID in case there a proteins with matches with more than one motif
for motif_n in ELM_keys:#for each ELM motif
    for ID in keys: #for every protein in the SP_ProteinList dict
        motif_m = 0 #counter to know the number of matches
        for m in re.finditer(ELM_models[motif_n][0], info_proteins[ID][2]):
            #print(motif_n)
            try:
                motif_i = m.group(1) #save the  instance
            except:
                print(m)
                continue
            motif_l = len(motif_i) #get the length
            motif_s = m.start() #get the start index

            motif_d = np.mean(proteins_dis[ID][motif_s:motif_s+motif_l]) #calculate the average disorder of the instance
            
            if motif_d > args.input_IUP_thrs:
                motif_c = "disorder"
            else:
                motif_c = "order"
            
            proteins_motifs[match_id] = [ID,motif_n, motif_m+1, motif_i ,motif_s,motif_d, motif_c] #save motif info: protein ID, name, match number, instance, start, average disorder
            match_id += 1 #change the ID for next match
            del motif_i, motif_l, motif_s, motif_d, motif_c #delete variable after storing their value
            motif_m += 1
            del m
        del motif_m

del match_id, motif_n


'''
Save Results in a text file
'''
print("Saving info.")

Results = list(proteins_motifs.keys())

os.chdir('.') #choose directory to save table
out_file = open(args.data_alias+"_MotifMatches_list.txt", "w") #name of the table
out_file.write('Protein_ID\tMotif_Name\tMatch_N\tMotif_Instance\tMotif_sSite\tMotif_Disorder\tDis_context\n' ) #header
for key in Results:
    out_file.write(f'{proteins_motifs[Results[key]][0]}\t{proteins_motifs[Results[key]][1]}\t{proteins_motifs[Results[key]][2]}\t{proteins_motifs[Results[key]][3]}\t{proteins_motifs[Results[key]][4]}\t{proteins_motifs[Results[key]][5]}\t{proteins_motifs[Results[key]][6]}\n' )
out_file.close()

del Results
