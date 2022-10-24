#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 13:57:14 2021

@author: JAVlvrd

Toxo Proteome and a list of motifs and motif hit sites

It will give as output a table with the conservation presence scores for each motif
"""

import argparse
import os
import pandas as pd
import re

'''
PATH OF PROTEOME
'''

parser=argparse.ArgumentParser()
parser.add_argument("input_ELMs",
                    help="path to a motif table file with names and REGEX") #receive file of signalP results
parser.add_argument("input_sites",
                    help="path to a motif match site table file with protein IDs and Motif sites") #receive file of signalP results
args = parser.parse_args()


#python MotifsInAlignments.py input_ELMs input_sites
#sample command
    #python ../Software/MotifMatches_InAlignments.py ../Data/elm_classes_200922.tsv ELM_Sep22_MotifMatches_sites.txt 

#outputs:
    #ALIAS_MotifMatches_presence.txt


strains = ['TGGT1',  'TGMAS', 'TGVAND',  'TGVEG' ]

species = ['BESB' , 'CSUI' ,  'HHA' ,'NCLIV' ]

out_dir = os.getcwd()

'''
Load ELM models

'''

ELMs = args.input_ELMs #Table with target ELM models
#ELMs = '/Users/JAVlvrd/Documents/Toxoplasma-2022/ToxoMotifs/Data/elm_classes_200922.tsv'
ELM_table = open(ELMs, 'r') #open up file
del ELMs

ELM_models = {}

for line in ELM_table:
    row = line.strip().split() #Separate the name and the regular expression
    Name = row[0] #select the Name of the ELM
    REGEX = row[1] #select the REGEX of the ELM
    Group = row[2] #select the slim group

    ELM_models[Name]=[REGEX, Group]
    del row
del ELM_table,line, Name, Group #delete used variables

ELM_keys = list(ELM_models.keys()) #Get ELM dict keys


'''
Load motif matches info

'''
MotifMatches_sites = args.input_sites
#MotifMatches_sites = '/Users/JAVlvrd/Documents/Toxoplasma-2022/ToxoMotifs/Results/ELM_Sep22_MotifMatches_sites.txt'
MotifMatches_sites_table = pd.read_table(MotifMatches_sites)
MotifMatches_sites_list = MotifMatches_sites_table.transpose().to_dict('list') #put the datafram in dictionary (of lists) format
del MotifMatches_sites_table

for key in  list(MotifMatches_sites_list.keys()):
    MotifMatches_sites_list[key][3] = list(map(int,MotifMatches_sites_list[key][3].split(",")))



'''
Read Alignments and get motif matches positions
'''

os.chdir('/Users/JAVlvrd/Documents/Toxoplasma-2021/Motif_Enrichments/Fastas/MassSpec_genes')

Aln_groups={}

for group in  list(MotifMatches_sites_list.keys()):
    
    ID=MotifMatches_sites_list[group][0] #TGME49_254720-t26_1-p1
    ID = re.sub("-t26_1-p1", "", ID)
    motif=MotifMatches_sites_list[group][1]


    aln_name = ID+'_groupIDs.txt.aln'
    try:
         aln_txt = open(aln_name, 'r')
    except:
        continue

    info_subjects = {} #dict (ID) will contain organism, sequence, seq length, disorder
    #slims_seqs = {}

    for line in aln_txt:
        if line.startswith('>'): #identify the lines with header with '>'
            header = line.strip().split('|') #Separate the different pieces of information and store them in a list

            subject_seqID = header[0] #select the ID of the Sequence
            subject_seqID = subject_seqID.replace('>','').replace(' ','') #take away extra characters and spaces

            info_subjects[subject_seqID]=['',0] #crate an entry in the info_proteins dic. using ID as key and storing gene name, aa length and empty sequence
            #adding an empty string to store amino acid sequence by extension

        else: #if the line is not a header then add it to the previous entry seq string
            aln_piece = line.strip() #take away extra formating from the line with aa's
            info_subjects[subject_seqID][0] = info_subjects[subject_seqID][0] + aln_piece #add currrent line to extend the protein sequence

            del aln_piece #delete
    del line, header, subject_seqID #delete used variables


    subject_keys = list(info_subjects.keys()) #Get ELM dict keys
    aln_group = {}

    for seq in subject_keys:    #getting real positions and alignment positions
        subject_seq=info_subjects[seq][0]
        subject_ref=[]
        counter=0
        for position  in subject_seq:
            if position == '-':
                subject_ref.append('-')
            else:
                subject_ref.append(counter)
                counter += 1
        del position

        slim_model=ELM_models[motif][0]
        slim_indexes = []


        #slim = re.findall(slim_model, subject_seq.replace("-",""))
        slim_index = re.finditer(slim_model, subject_seq.replace("-",""))
        for m in slim_index:
            #print(motif_n)
            #motif_i = m.group(1) #save the  instance
            #motif_l = len(motif_i) #get the length
            motif_s = subject_ref.index(m.start()) #get the start index
            slim_indexes.append(motif_s)
            del motif_s #,motif_i,motif_l

        if seq.split("_")[0] == "TGME49" :
            aln_group[seq+".ref"] = list(subject_ref[x] for x in slim_indexes)
            aln_group[seq] = slim_indexes
        else:
            aln_group[seq] = slim_indexes

    del seq,subject_seq, subject_ref, counter, slim_model,slim_indexes,slim_index, m
    Aln_groups[ID+"|"+motif] = aln_group
del group, motif, aln_name, aln_txt,subject_keys, info_subjects, aln_group


# for group in  list(MotifMatches_sites_list.keys()):
#     ID = MotifMatches_sites_list[group][0] #TGME49_254720-t26_1-p1
#     ID = re.sub("-t26_1-p1", "", ID)
#     motif=MotifMatches_sites_list[group][1]
#     key = ID+"|"+motif

#     if motif == "PTAP":
#         print(MotifMatches_sites_list[group])
#         for seq in list(Aln_groups[key].keys()):
#             print(group,key, seq, Aln_groups[key][seq])
#         print("\n")
# del group


DisOrgComp_cons = {}
motif_window = 15

for group in  list(MotifMatches_sites_list.keys()):
    ID = MotifMatches_sites_list[group][0] #TGME49_254720-t26_1-p1
    ID = re.sub("-t26_1-p1", "", ID)
    motif=MotifMatches_sites_list[group][1]
    key = ID+"|"+motif


    #if motif == "PTAP":
    true_motifs = MotifMatches_sites_list[group][3]

    try:
        all_motifs = Aln_groups[key][ID+'-t26_1-p1.ref'] #when there's no matches of a motif in a certain alignment e.g
    except:
        continue

    true_index =[]
    for site in true_motifs:
        true_index.append(all_motifs.index(site))
    del site

    guide_index = list(Aln_groups[key][ID+'-t26_1-p1'][x] for x in true_index)
    #print(MotifMatches_sites_list[group], true_index,guide_index, "\n")


    for ind in guide_index:
        seq_n = 0
        strains_n = 0
        species_n = 0

        ind_group = {k: value for k, value in Aln_groups[key].items() if ind in value}
        #ind_seq = MotifMatches_sites_list[group][2][guide_index.index(ind)]
        ind_seq = MotifMatches_sites_list[group][2]

        for seq in list(ind_group.keys()):
            interval = range(ind-motif_window, ind+motif_window)
            #if ind in ind_group[seq]:
            if bool(set(interval) & set(ind_group[seq])):  #here should be a more flexible criterium of motif position
                org = seq.split("_")[0]
                if org in strains:
                    strains_n += 1
                elif org in species:
                    species_n += 1

                seq_n += 1

        cons_st = strains_n/4
        cons_sp = species_n/4
        cons_org = seq_n/9
        #print(key, ind,cons_org,cons_st,cons_sp)
        DisOrgComp_cons[key+"|"+str(ind_seq)] = [cons_org,cons_st,cons_sp] #correct index to position
        del seq, org, interval
        del ind, seq_n,strains_n, species_n,ind_group,ind_seq,cons_org,cons_st,cons_sp
del true_motifs, all_motifs, true_index,guide_index,group
del ID, motif, key



#species_t = motif_sp_n / species_n
#species_r = motif_sp_n / species_n


'''
Save presences scores in a text file
'''

print('Writing output file')

os.chdir(out_dir) #choose directory to save table
out_name = re.sub("_sites.txt", "", MotifMatches_sites)
out_file = open(out_name+"_presence.txt", "w") #name of the table
out_file.write('key\t seq_id\t motif\t motif_site\tpresence_org\t presence_str\t presence_spc\n')#header
for key in list(DisOrgComp_cons.keys()):
    info = key.split('|')
    out_file.write(f'{key}\t{info[0]}\t{info[1]}\t{info[2]}\t{DisOrgComp_cons[key][0]}\t{DisOrgComp_cons[key][1]}\t{DisOrgComp_cons[key][2]}\n' )
out_file.close()
del key, info, out_file
