#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 13:57:14 2021

@author: JAVlvrd

Toxo Proteome and a list of motifs and motif hit sites

It will give as output a table with the conservation presence scores for each motif
"""



'''
PATH OF PROTEOME
'''

import argparse
import os
import pandas as pd
import re

parser=argparse.ArgumentParser()
parser.add_argument("input_list",
                    help="path to the motif match list") #receive file of signalP results
parser.add_argument("input_PTMs",
                    help="path to a PTM table file with protein ID, aa and residue number") #receive file of signalP results
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
MotifMatches_table = pd.read_table(MotifMatches)
MotifMatches_list = MotifMatches_table.transpose().to_dict('list') #put the datafram in dictionary (of lists) format
del MotifMatches, MotifMatches_table


for key in  list(MotifMatches_list.keys()):    
    seq_id = ID = re.sub("-t26_1-p1", "", MotifMatches_list[key][0])
    MotifMatches_list[key].append(seq_id)
del seq_id, key


'''
Load PTM info

'''

PTMs_file = args.input_PTMs
#PTMs_file = '/Users/JAVlvrd/Documents/Toxoplasma-2021/Motif_Enrichments/ME49_Phosphosites_130622.txt_byModification.tab'
PTMs_table = pd.read_table(PTMs_file)
PTMs_table.columns = PTMs_table.columns.str.replace(" ","") 
del PTMs_file



'''
Count and add PTM info to Motif matches list

'''

MotifMatches_PTMs = {}

ptm_n = 0
for ptm in list(MotifMatches_list.keys()):
    seq_id = MotifMatches_list[ptm][7]
    motif_c = MotifMatches_list[ptm][1]
    m_start = MotifMatches_list[ptm][4] + 1 #to correct for indexing
    m_ins = MotifMatches_list[ptm][2]
    m_end = m_start + len(MotifMatches_list[ptm][3])
    
    mods = PTMs_table[ (PTMs_table["seq_id"] == seq_id) & (PTMs_table["residue"] >= m_start) & (PTMs_table["residue"] <= m_end) ]
    
    if len(mods) > 0:
        mods_n = len(mods)
        mods_res = mods['residue'].apply(str).to_list()
        mods_res = mods.aminoacid + mods_res
        mods_res = ",".join(mods_res.to_list())
        
        MotifMatches_PTMs[ptm_n] = [seq_id,motif_c,m_ins,m_start,mods_n,mods_res]
        ptm_n += 1
        
        del mods_n, mods_res
       
    del mods, ptm
del seq_id, motif_c, m_ins, m_start, m_end, ptm_n         


'''
Save PTM-motif list in a text file
'''

print('Writing output file')

os.chdir(out_dir) #choose directory to save table
out_name = re.sub("_list.txt", "", args.input_list)
out_file = open(out_name+"_mods.txt", "w") #name of the table
out_file.write('key\tmodsNum\tmodsSites\n')#header
for key in list(MotifMatches_PTMs.keys()):
    motif_id = MotifMatches_PTMs[key][0] + '|' + MotifMatches_PTMs[key][1] + '|' + str(MotifMatches_PTMs[key][2])
    out_file.write(f'{motif_id}\t{MotifMatches_PTMs[key][4]}\t{MotifMatches_PTMs[key][5]}\n')
out_file.close()
del key, motif_id, out_file
del out_name

del MotifMatches_PTMs
del MotifMatches_list,PTMs_table
