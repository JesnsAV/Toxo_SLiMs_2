#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 2022

@author: JAVlvrd


"""



'''
INPUT & ARGUMENTS
'''

import argparse
import os
import pandas as pd
import re

parser=argparse.ArgumentParser()
parser.add_argument("input_list",
                    help="path to the motif match list") #receive file of motif
parser.add_argument("input_Domains",
                    help="path to a Domain table from Uniprot with protein ID, domain and residue numbers") #receive file of signalP results
parser.add_argument("input_mappings",
                    help="path to an ID mapping table between Uniprot and ToxoDB") #receive file of signalP results
args = parser.parse_args()


#python MotifsInAlignments.py input_ELMs input_sites
#sample command
    #

#outputs:
    #


out_dir = os.getcwd()

'''
Load Motif match lists 

'''

MotifMatches = args.input_list
#MotifMatches = '/Users/JAVlvrd/Documents/Toxoplasma-2022/ToxoMotifs/Results/ELM_Sep22_MotifMatches_list.txt'
MotifMatches = '/Users/JAVlvrd/Documents/Toxoplasma-2022/ToxoMotifs/Results/ELM_Dec22/ELM_Dec22_MotifMatches_list.txt'
MotifMatches_table = pd.read_table(MotifMatches)
MotifMatches_table['seq_id'] = MotifMatches_table['Protein_ID'].str.replace("-t26_1-p1","")
MotifMatches_table['Motif_eSite'] = MotifMatches_table['Motif_sSite'] + MotifMatches_table['Motif_Instance'].apply(len)
#MotifMatches_table['key'] = MotifMatches_table['seq_id'] + "|" + MotifMatches_table['Motif_Name'] + "|"  + MotifMatches_table['Match_N'].to_string()
MotifMatches_table['key'] = MotifMatches_table['seq_id'] + "|" + MotifMatches_table['Motif_Name'] + "|"  + MotifMatches_table['Match_N'].astype(str)
#MotifMatches_list = MotifMatches_table.transpose().to_dict('list') #put the datafram in dictionary (of lists) format
del MotifMatches


# for key in  list(MotifMatches_list.keys()):    
#     seq_id = ID = re.sub("-t26_1-p1", "", MotifMatches_list[key][0])
#     MotifMatches_list[key].append(seq_id)
# del seq_id, key


'''
Load Domain info

'''

DOMs_file = args.input_Domains
#DOMs_file = '/Users/JAVlvrd/Documents/Toxoplasma-2022/ToxoMotifs/Data/uniprot-dom-2022.11.18-12.00.54.21.tsv'
DOMs_table = pd.read_table(DOMs_file)
DOMs_table.columns = DOMs_table.columns.str.replace(" ","") 
del DOMs_file


'''
Load mapping info

'''

MAPs_file = args.input_mappings
#MAPs_file = '/Users/JAVlvrd/Documents/Toxoplasma-2022/ToxoMotifs/Data/uniprot-mapp-2022.07.29-13.13.55.06.tsv'
MAPs_table = pd.read_table(MAPs_file)
MAPs_table.columns = MAPs_table.columns.str.replace(" ","") 
del MAPs_file



'''
Add ToxoDB ID to Domain info

'''
print('Adding domain info')
DOMs_a_table = pd.merge(DOMs_table, MAPs_table, on="From")
DOMs_a_table.set_index('From')
DOMs_list = DOMs_a_table.transpose().to_dict('list')
del DOMs_a_table, DOMs_table, MAPs_table

DOMs = {}
for key in list(DOMs_list.keys()):
    ID = DOMs_list[key][3].replace("ToxoDB:","") 
    if isinstance(DOMs_list[key][2], str):
        DOMs[ID] = [DOMs_list[key][0],DOMs_list[key][2]]
del key, ID
del DOMs_list

MotifMatches_table['doms_num'] = 0
MotifMatches_table['doms_name'] = 'NA'

for dom in list(DOMs.keys()):
    #prot = 'TGME49_227030'
    prot_motifs = MotifMatches_table.loc[MotifMatches_table['seq_id'] == dom]
    if len(prot_motifs) > 0:
        doms = DOMs[dom][1].split("DOMAIN")
        for i in range(len(doms)-1):
            info_doms = doms[i+1].split(";")
            pos_doms = info_doms[0].replace(' ','').split("..")
            name_doms = info_doms[1].replace(' /note=','').replace('\"','')
            motif_doms = prot_motifs.loc[ (prot_motifs['Motif_sSite'] > int(pos_doms[0])) & (prot_motifs['Motif_eSite'] < int(pos_doms[1])) ]
            MotifMatches_table.iloc[motif_doms.index,[10]] += 1
            MotifMatches_table.iloc[motif_doms.index,[11]] = MotifMatches_table.iloc[motif_doms.index,[11]] + name_doms
#     del i
 
# del prot_motifs, doms, info_doms, pos_doms, name_doms,motif_doms
# del dom

# del DOMs

'''
Save PTM-motif list in a text file
'''

print('Writing output file')

os.chdir(out_dir) #choose directory to save table
out_name = re.sub("_list.txt", "", args.input_list)
out_file = out_name+"_doms.txt" #name of the table
MotifMatches_DOMs = MotifMatches_table[["key","doms_num","doms_name"]]
MotifMatches_DOMs.to_csv('ELM_Dec22_MotifMatches_doms.txt', sep="\t")

del MotifMatches_table,MotifMatches_DOMs







