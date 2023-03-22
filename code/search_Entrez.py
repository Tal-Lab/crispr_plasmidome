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
import metadata

#Restrict request to only ask for HTTP/1.0
http.client.HTTPConnection._http_vsn = 10
http.client.HTTPConnection._http_vsn_str = 'HTTP/1.0'

#Read BLASTn with plasmids metadata taxonomy

BLASTn_with_metadata = metadata.plasmid_blast_copy

#Creating new data frame with unique spacer id

BLASTn_with_metadata2 = BLASTn_with_metadata[["sseqid"]]
BLASTn_with_metadata2 = BLASTn_with_metadata2.drop_duplicates()

Entrez.email = "androsiu@post.bgu.ac.il"  # Tell NCBI who you are

#Pull specise from NCBI using Entrez.
#In order to pull the rest of the taxonomy, change "organism" to "taxonomy".
def search_Entrez (x):
    try:
        handle = Entrez.efetch(db="nucleotide", id= x, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        return [record.annotations["taxonomy"], record.annotations["organism"]]
    except ValueError:
        return "Error"

#Apply function on id's column from the file
BLASTn_with_metadata2 [["taxonomy", "species"]]= BLASTn_with_metadata2 ["sseqid"].apply(search_Entrez, axis = 1, result_type='expand')

#Merge BLASTn with plasmids metadata taxonomy and spacers taxonomy
BLASTn_with_metadata = BLASTn_with_metadata.merge(BLASTn_with_metadata2, on='sseqid')

BLASTn_with_metadata.to_csv("BLAST_DataBase.csv")
