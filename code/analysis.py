# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 14:40:00 2022

@author: sivan
"""


import pandas as pd
import glob
import os

#Read BLASTn file
final_BLASTp = pd.read_csv (r"C:\sivan\לימודים\שנה ג\פרויקט מחקר בפייתון/BLASTp_DataBase3.csv", sep= ",",  header = 0)

#Count number of occurrences of each values in "match" column 
sum_df  = final_BLASTp['match'].value_counts()

#Isolate match by species level
species_unique = final_BLASTp.loc[final_BLASTp['match'] == 'species']
#isolata unique ids with species match in a new dataframe
unique_species_match = pd.DataFrame(species_unique['qseqid'].unique())
#Add column "match" with species value to all unique id
unique_species_match['match'] = 'species'

#Isolate match by genus level
genuse_unique = final_BLASTp.loc[final_BLASTp['match'] == 'genus']
#isolata unique ids with genus match in a new dataframe
unique_match = pd.DataFrame(genuse_unique['qseqid'].unique())
#Add column "match" with genus value to all unique id
unique_match['match'] = 'genus'
#Filter genus results by removing results that also have species match   
unique_match_genus = unique_match [~unique_match[0].isin(unique_species_match[0])]

#Isolate match by famaliy level
family_unique = final_BLASTp.loc[final_BLASTp['match'] == 'family']
#Isolata unique ids with family match in a new dataframe
unique_family_match = pd.DataFrame(family_unique['qseqid'].unique())
#Add column "match" with family value to all unique id
unique_family_match["match"] = 'family'
#Filter family results by removing results that also have species match   
unique_match_family_species = unique_family_match [~unique_family_match[0].isin(unique_species_match[0])]
#Filter family results by removing results that also have genus match   
unique_match_family = unique_match_family_species [~unique_match_family_species[0].isin(unique_match[0])]

#Isolate match by order level
order_unique = final_BLASTp.loc[final_BLASTp['match'] == 'order']
#Isolata unique ids with order match in a new dataframe
unique_order_match = pd.DataFrame(order_unique['qseqid'].unique())
#Add column "match" with order value to all unique id
unique_order_match["match"] = 'order'
#Filter order results by removing results that also have species match   
unique_match_order_species = unique_order_match [~unique_order_match[0].isin(unique_species_match[0])]
#Filter order results by removing results that also have genus match   
unique_match_order_genus = unique_match_order_species [~unique_match_order_species[0].isin(unique_match[0])]
#Filter order results by removing results that also have family match   
unique_match_order_family = unique_match_order_genus[~unique_match_order_genus[0].isin(unique_match_family[0])]

#Isolate match by class level
group_class = final_BLASTp.loc[final_BLASTp['match'] == 'class']
#Isolata unique ids with calss match in a new dataframe
unique_class_match = pd.DataFrame(group_class['qseqid'].unique())
#Add column "match" with class value to all unique id
unique_class_match["match"] = 'class'
#Filter class results by removing results that also have species match   
unique_match_class_species = unique_class_match [~unique_class_match[0].isin(unique_species_match[0])]
#Filter class results by removing results that also have genus match 
unique_match_class_genus = unique_match_class_species[~unique_match_class_species[0].isin(unique_match[0])]
#Filter class results by removing results that also have family match   
unique_match_class_family = unique_match_class_genus[~unique_match_class_genus[0].isin(unique_match_family[0])]
#Filter class results by removing results that also have order match   
unique_match_class_order = unique_match_class_family[~unique_match_class_family[0].isin(unique_match_order_family[0])]

#Isolate match by phylum level
group_phylum = final_BLASTp.loc[final_BLASTp['match'] == 'phylum']
#Isolata unique ids with phylum match in a new dataframe
unique_phylum_match = pd.DataFrame(group_phylum['qseqid'].unique()) 
#Add column "match" with phylum value to all unique id
unique_phylum_match["match"] = 'phylum'
#Filter phylum results by removing results that also have species match   
unique_match_phylum_species = unique_phylum_match [~unique_phylum_match[0].isin(unique_species_match[0])]
#Filter phylum results by removing results that also have genus match 
unique_match_phylum_genus = unique_match_phylum_species[~unique_match_phylum_species[0].isin(unique_match[0])]
#Filter phylum results by removing results that also have family match  
unique_match_phylum_family = unique_match_phylum_genus[~unique_match_phylum_genus[0].isin(unique_match_family[0])]
#Filter phylum results by removing results that also have order match  
unique_match_phylum_order = unique_match_phylum_family[~unique_match_phylum_family[0].isin(unique_match_order_family[0])]
#Filter phylum results by removing results that also have class match
unique_match_phylum_class = unique_match_phylum_order[~unique_match_phylum_order[0].isin(unique_match_class_order[0])]

#Isolate match by superkingdom level
group_kingdom = final_BLASTp.loc[final_BLASTp['match'] == 'superkingdom']
#Isolata unique ids with superkingdom match in a new dataframe
unique_kingdom_match = pd.DataFrame(group_kingdom['qseqid'].unique())
#Add column "match" with superkingdom value to all unique id 
unique_kingdom_match["match"] = 'superkingdom'
#Filter phylum results by removing results that also have species match   
unique_match_kingdom_species = unique_kingdom_match [~unique_kingdom_match[0].isin(unique_species_match[0])]
#Filter superkingdom results by removing results that also have genus match 
unique_match_kingdom_genus = unique_match_kingdom_species[~unique_match_kingdom_species[0].isin(unique_match[0])]
#Filter superkingdom results by removing results that also have family match 
unique_match_kingdom_family = unique_match_kingdom_genus[~unique_match_kingdom_genus[0].isin(unique_match_family[0])]
#Filter superkingdom results by removing results that also have order match 
unique_match_kingdom_order = unique_match_kingdom_family[~unique_match_kingdom_family[0].isin(unique_match_order_family[0])]
#Filter superkingdom results by removing results that also have class match
unique_match_kingdom_class = unique_match_kingdom_order[~unique_match_kingdom_order[0].isin(unique_match_class_order[0])]
#Filter superkingdom results by removing results that also have phylum match
unique_match_kingdom_phylum = unique_match_kingdom_class[~unique_match_kingdom_class[0].isin(unique_match_phylum_class[0])]

#Merge superkindom data with analysis result in ordr to find reasons for this findings
unique_match_kingdom_phylum = unique_match_kingdom_phylum.merge(final_BLASTp, left_on=0, right_on='qseqid')
#Searching for known host defind as 'uncultured bacterium'
uncultured_bacterium = unique_match_kingdom_phylum.loc[((unique_match_kingdom_phylum['plasmid species'] == 'uncultured bacterium') | (unique_match_kingdom_phylum['plasmid species'] == 'uncultured bacterium HHV216') | (unique_match_kingdom_phylum['plasmid species'] == 'uncultured bacterium pMCBF6')) ]
#Finding unique uncultured bacterium
unique_uncultured_bacterium = uncultured_bacterium['qseqid'].unique()

