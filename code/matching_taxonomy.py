# -*- coding: utf-8 -*-
"""
Created on Wed May  4 11:33:51 2022

@author: sivan
"""


from Bio import SeqIO
from Bio import Entrez
import pandas as pd
import sys
import http.client
import numpy as np


final_BLASTp = pd.read_csv ("/gpfs0/tals/projects/Analysis/sivan_project/BLASTp_DataBase2.csv", sep= ",",  header = 0)

#Turn the lists in "spacer host taxonomy" column in to type list
#First removing all unnecessary elements
final_BLASTp["spacer host taxonomy"] = final_BLASTp["spacer host taxonomy"].str.replace("'","")

final_BLASTp["spacer host taxonomy"] = final_BLASTp["spacer host taxonomy"].str.replace("[","")

final_BLASTp["spacer host taxonomy"] = final_BLASTp["spacer host taxonomy"].str.replace("]","")

final_BLASTp["spacer host taxonomy"] = final_BLASTp["spacer host taxonomy"].str.replace(" ","")

#Making a list from values in a cell
def to_list (x):
    return list(x.split( ","))

#Apply for all values in "spacer host taxonomy" column
final_BLASTp["spacer host taxonomy"] = final_BLASTp["spacer host taxonomy"].apply(to_list)




#Finding a match between the taxonomy of the spacer and the taxonomy of the plasmid in each alignment result
for i in range(1,9):
    for index, row in final_BLASTp.iterrows():
        if (pd.isna(final_BLASTp['match'][index]) == True):
            if (final_BLASTp['spacer host species'][index] == final_BLASTp['plasmid species'][index]) == True:
                final_BLASTp['match'][index] = 'species'
                continue
            elif (final_BLASTp['spacer host taxonomy'][index][-i] == final_BLASTp['plasmid genus'][index]) == True:
                final_BLASTp['match'][index] = 'genus' 
                continue
            elif (final_BLASTp['spacer host taxonomy'][index][-i] == final_BLASTp['plasmid family'][index]) == True:
                final_BLASTp['match'][index] = 'family'
                continue 
            elif (final_BLASTp['spacer host taxonomy'][index][-i] == final_BLASTp['plasmid order'][index]) == True:
                final_BLASTp['match'][index] = 'order'
                continue 
            elif (final_BLASTp['spacer host taxonomy'][index][-i] == final_BLASTp['plasmid class'][index]) == True:
                final_BLASTp['match'][index] = 'class' 
                continue
            elif (final_BLASTp['spacer host taxonomy'][index][-i] == final_BLASTp['plasmid phylum'][index]) == True:
                final_BLASTp['match'][index] = 'phylum' 
                continue
            elif (final_BLASTp['spacer host taxonomy'][index][-i] == final_BLASTp['plasmid superkingdom'][index]) == True:
                final_BLASTp['match'][index] = 'superkingdom'
                continue
        else:
            continue
            
            
print("done")
   
           
final_BLASTp.to_csv("BLASTp_DataBase3.csv")
