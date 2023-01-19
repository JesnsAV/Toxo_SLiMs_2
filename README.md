# Toxo_SLiMs_2
Finding SLiMs in Toxoplasma gondii proteome

Complete table with raw matches
ELM_Dec22_MotifMatches_complete.tsv.gz

Complete table with filtered matches
ELM_Dec22_MotifMatches_filtered.tsv.gz


#In order to produced the previous tables

###FIRST MATCH DETERMINATION
##INPUTS
#ELM motif classes table format:
#CLV_C14_Caspase3-7 ([DSTE][^P][^DEWHFYC]D[GSAN]) 0.00309374

#ToxoDB proteome in fasta format:
#>TGME49_287280-t26_1-p1
#MMHLIQKKCPGFPPGFQLPCRLKARRGRLFRHESCTMLFSVALCLTALASFVPFECSTRR…
#Command for motif match search:
Python MotifMatches_Dis.py 'ALIAS' TgondiiME49.fasta Elm_classes.tsv 0.4

##OUPUT
#Motif matches result table format:
#Protein_ID Motif_Name Match_N Motif_Instance Motif_sSite Motif_Disorder Dis_context
#TGME49_287280-t26_1-p1 CLV_C14_Caspase3-7 1 TERDG 85 0.634 disorder

###Reformatting motif match sites:
Rscript MotifMatches_Sites.R ALIAS_MotifMatches_list.txt

###Integrating alignments info
##INPUTS
#Motif match site table format:
#Protein_ID Motif_Name Match_N Motif_sites
#TGME49_287280-t26_1-p1 CLV_C14_Caspase3-7 1 85 109

##Command to find motif matches in alignments:
Python MotifMatches_InAlignments.py Elm_classes.tsv ALIAS_MotifMatches_sites.txt

##OUPUT
#Motif match presence table format:
#key seq_id motif motif_site presence_org presence_str presence_spc
#TGME49_292920|CLV_C14_Caspase3-7|1 TGME49_292920 CLV_C14_Caspase3-7 1 0.7777 1.0 0.5

###Integrating AlphaFold and DSSP information
##Command to map pLDDT and accessibility values to motif matches:
Python MotifMatches_pLDDT.py ALIAS_MotifMatches_list.txt TgondiiME49_AF_plddt_values.txt TgondiiME49_AF_dssp_values.txt

##OUTPUT
#Motif match pLDDT and accessibility values table format:
#Key pLDDT Accessibility
#TGME49_287280|CLV_C14_Caspase3-7|1 31.3800 0.7203

###Integrating Phosphosite information
##Command to phosphosite to motif matches:
Python MotifMatches_PTMs.py ALIAS_MotifMatches_list.txt TgondiiME49_Phosphosites.tab

##OUTPUT
#Motif match phosphosite table format:
#Key modsNum modsSites
#TGME49_293300|CLV_C14_Caspase3-7|1 2 S96, S99

###Integrating Domain mappings
##Command to add domain mapping to motif matches:
Python MotifMatches_Domains.py ALIAS_MotifMatches_list.txt TgondiiME49_DomainMappings.tsv TgondiiME49_DBMappings.tsv

##OUTPUT
#Motif match domain mapping table format:
#key doms_num doms_name
#TGME49_305460|DOC_PP2A_B56_1|1 0 NA
#TGME49_305460|DOC_PP2A_B56_1|2 1 NAPeptidase_M24

###Integrating all information together
#Command to combine all the motif matches information:
Rscript MotifMatches_Enrichment.R ALIAS

##OUTPUT
#Protein_ID	Motif_Name	Match_N	Motif_Instance	Motif_sSite	Motif_Disorder	Dis_context	seq_id	key	Motif_Type	presence_org	presence_str	presence_spc	ExpEvd	Product_Description	Final_Probability_TAGM_MAP	Predicted_Location_TAGM_MAP	TM_Domains	LocEvd	Organelle	modsNum	modsSites	...1	doms_num	doms_name	mean_pLDDT_a	mean_acc_a	mean_pLDDT_b	mean_acc_b
#TGME49_287280-t26_1-p1	CLV_C14_Caspase3-7	1	TERDG	85	0.6344230000000001	disorder	TGME49_287280	TGME49_287280|CLV_C14_Caspase3-7|1	CLV	NA	NA	NA	no	NA	NA	NA	NA	no	NA	NA	NA	0	0	NA	31.380000000000003	0.7203526698508012	NA	NA

###To Produce Filtered results
#MotifMatches_Tables.Rmd
