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
import os, re
from pathlib import Path

#Restrict request to only ask for HTTP/1.0
http.client.HTTPConnection._http_vsn = 10
http.client.HTTPConnection._http_vsn_str = 'HTTP/1.0'

# Cluster
path = r"/gpfs0/tals/projects/Analysis/Lucy_plasmidome/Plasmidome/CRISPR"

# working directories
visuals = f"{path}/visualisations"
tables = f"{path}/data_calculations"
resource = f"{path}/res"
Path(visuals).mkdir(parents=True, exist_ok=True)
Path(tables).mkdir(parents=True, exist_ok=True)

# working files
unique_spacers = f"{resource}/unique_spacers.csv"

#Read BLASTn with plasmids metadata taxonomy

BLASTn_with_metadata = pd.read_csv(unique_spacers, sep= ",",  header = 0, index_col = 0)
BLASTn_with_metadata = BLASTn_with_metadata.drop_duplicates()

Entrez.email = "maane@post.bgu.ac.il"  # Tell NCBI who you are

#Pull specise from NCBI using Entrez.
#In order to pull the rest of the taxonomy, change "organism" to "taxonomy".
def search_Entrez (x):
    try:
        handle = Entrez.efetch(db="nucleotide", id= x, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        return record.description
    except ValueError:
        return "Error"

#Apply function on id's column from the file
BLASTn_with_metadata ["description"]= BLASTn_with_metadata["unique"].apply(search_Entrez)

def search_plasmid_or_not (x):
    try:
        find = 'plasmid'
        if find in x:
            return find
        else:
            return "not plasmid"
    except ValueError:
        return "Error"

BLASTn_with_metadata ["pplasmid_or_not"]= BLASTn_with_metadata ["description"].apply(search_plasmid_or_not)
BLASTn_with_metadata.to_csv(f"{tables}/plasmid_or_not.csv")