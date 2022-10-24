#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 16:02:52 2022

@author: JAVlvrd
"""


'''
PATH OF PROTEOME
'''

import argparse
import re
import glob
import numpy as np
import os

parser=argparse.ArgumentParser()
parser.add_argument("input_list",
                    help="Table with the list of motif matches ") #receive file of signalP results
parser.add_argument("AF_values",
                    help="Table with pLDDT alphafold values") #receive proteome path 
parser.add_argument("DSSP_values",
                    help="Table with accessibility DSSP  values ") #receive proteome path 
args = parser.parse_args()

#sample command
    #python ../Software/MotifMatches_pLDDT.py ELM_test_MotifMatches_list.txt ../Data/TgondiiME49_AF_plddt_values.txt ../Data/TgondiiME49_AF_dssp_values.txt 

#outputs:
    #ELM_test_MotifMatches_AF_extension.txt


out_dir = os.getcwd()

'''
Load pLDDT and DSSP accesibility values
'''
import sys
import pandas as pd

plddt_file = args.AF_values
#plddt_file = '/Users/JAVlvrd/Documents/Toxoplasma-2022/ToxoMotifs/Data/TgondiiME49_AF_plddt_values.txt'
plddt_table = pd.read_table(plddt_file)
plddt_list_raw = plddt_table.set_index('AF_ID').transpose().to_dict('list') #put the datafram in dictionary (of lists) format
plddt_list = {}
for item in plddt_list_raw:
    temp = plddt_list_raw[item][0]
    temp = re.sub("\[|\]", "", temp)
    temp = temp.split(",")
    temp = [ float(n) for n in temp ]
    plddt_list[item] = temp
del temp, plddt_list_raw
del plddt_file, plddt_table



access_file = args.DSSP_values
#access_file = '/Users/JAVlvrd/Documents/Toxoplasma-2022/ToxoMotifs/Data/TgondiiME49_AF_dssp_values.txt'
access_table = pd.read_table(access_file)
access_list_raw = access_table.set_index('AF_ID').transpose().to_dict('list') #put the datafram in dictionary (of lists) format
access_list = {}
for item in access_list_raw:
    temp = access_list_raw[item][0]
    temp = re.sub("\[|\]", "", temp)
    temp = temp.split(",")
    temp = [ float(n) for n in temp ]
    access_list[item] = temp
del temp, access_list_raw
del access_file, access_table

 

'''
Get ToxoDB ids
'''
ToxoDB_ID = {}

ToxoIDs = "/Users/JAVlvrd/Documents/Toxoplasma-2022/AlphaFold/AlphaFoldDB_Structures/uniprot-compressed_true_download_true_format_tsv-2022.07.29-13.13.55.06.tsv"
ToxoDB_ID_table = open(ToxoIDs, 'r') #open up file
del ToxoIDs


for line in ToxoDB_ID_table:
    row = line.strip().split() #Separate the name and the regular expression
    UP_id = row[0] #select the Name of the ELM
    Toxo_id = re.sub('ToxoDB:','',row[1]) #select the REGEX of the ELM
    ToxoDB_ID[Toxo_id]=[UP_id] 
    del row
del ToxoDB_ID_table,line, UP_id, Toxo_id #delete used variables
del ToxoDB_ID['To']     


'''
Get Motif Table
'''
Motif_file = args.input_list
#Motif_file = "/Users/JAVlvrd/Documents/Toxoplasma-2022/ToxoMotifs/Results/ELM_Sep22_MotifMatches_list.txt"
Motif_table = pd.read_table(Motif_file)
Motif_dict = Motif_table.transpose().to_dict('list') #put the datafram in dictionary (of lists) format
del Motif_table


'''
Add pLDDT values to motifs
'''

for key in list(Motif_dict.keys()):  
    seq_id = Motif_dict[key][0]
    seq_id = re.sub("-t26_1-p1","",seq_id)
    
    try:
        pdb_id = ToxoDB_ID[seq_id][0]
        pdb_id = 'AF-'+pdb_id+'-F1-model_v3.pdb'
        dssp_id = pdb_id+'.dssp'
        plddt_list[pdb_id]
        access_list[dssp_id]
    except:
        Motif_dict[key].append('na')
        Motif_dict[key].append('na')
        continue
    
    motif_st = Motif_dict[key][4]
    motif_ln = len(Motif_dict[key][3])
    motif_plddt_val = plddt_list[pdb_id][motif_st:motif_st+motif_ln]
    motif_acc_val = access_list[dssp_id][motif_st:motif_st+motif_ln]
          
    motif_plddt_scr = np.mean(motif_plddt_val)
    motif_acc_scr = np.mean(motif_acc_val)
    #print(key,seq_id,motif_score)
    #print(key)
    Motif_dict[key].append(motif_plddt_scr)
    Motif_dict[key].append(motif_acc_scr)
    
del key, seq_id,pdb_id,dssp_id
del access_list,plddt_list
del motif_st, motif_ln,motif_plddt_val,motif_acc_val, motif_plddt_scr, motif_acc_scr

'''
Save extended Motif Table
'''

os.chdir(out_dir) #choose directory to save table
#out_file = open(args.input_proteome+"_Disorder.txt", "w") #name of the table
out_name = re.sub("_list.txt", "", args.input_list)
out_file = open(out_name+"_AF_extension.txt", "w") #name of the table
for key in list(Motif_dict.keys()):
    seq_id = re.sub("-t26_1-p1","",Motif_dict[key][0])
    motif_id = seq_id + '|' + Motif_dict[key][1] + '|' + str(Motif_dict[key][2])
    out_file.write(f'{motif_id}\t{Motif_dict[key][-2]}\t{Motif_dict[key][-1]}\n' )
out_file.close()
del seq_id,motif_id
del key, out_file







