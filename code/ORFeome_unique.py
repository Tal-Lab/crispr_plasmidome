# -*- coding: utf-8 -*-
"""
Created on 17/11/2022 14:33

Author: Lucy
"""
'''
Retrieving ORFeome for unique plasmids with hits to spacers
'''

### importing modules
import pandas as pd
import numpy as np
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os, re
from pathlib import Path

#pd.set_option('display.max_colwidth', None)
### paths
email = "androsiu@post.bgu.ac.il"  # Tell NCBI who you are
# uncomment relevant path to OS
# Windows
path = r"C:\Users\Lucy\iCloudDrive\Documents\bengurion\Project students\Sivan_project"
# macOS
#path = r"/Users/lucyandrosiuk/Documents/bengurion/Project students/Sivan_project"
# Cluster
path = r"/gpfs0/tals/projects/Analysis/Lucy_plasmidome/Plasmidome/CRISPR"

# working directories
visuals = f"{path}/visualisations"
tables = f"{path}/data_calculations"
Path(visuals).mkdir(parents=True, exist_ok=True)

# working files
unique_plasmids = f"{path}/unique_plasmids.csv"

def search_Entrez (x):
    '''getting taxonomy information for each spacer-barer id'''
    Entrez.email = email  # Tell NCBI who you are
    Entrez.sleep_between_tries = 5 # Tell NCBI delay, in seconds, before retrying a request
    try:
        handle = Entrez.efetch(db="nuccore", id= x, rettype="fasta_cds_aa", retmode="text")
        records = SeqIO.parse(handle, "fasta")
        ids = []
        p_names = []
        seqs = []
        for record in records:
            if record.description.__contains__('[pseudo=true]'):
                pass
            else:
                id = record.id[4:]
                ids.append(id)
                description = record.description
                parts = description.split()
                p = parts[2:-3]
                p = " ".join(p)
                p_name = p[9:-1]
                p_names.append(p_name)
                sequence = str(record.seq)
                seqs.append(sequence)
        df = pd.DataFrame({'ID': ids, 'P_Name': p_names, 'Sequence': seqs})
        #print(record.annotations["structured_comment"])
        handle.close()
        return df
    except ValueError:
        return "Error"

"""
def search_Entrez (x):
    '''getting taxonomy information for each spacer-barer id'''
    Entrez.email = email  # Tell NCBI who you are
    Entrez.sleep_between_tries = 5 # Tell NCBI delay, in seconds, before retrying a request
    print("Entered search_Entrez function with %s" % x)
    handle = Entrez.efetch(db="nuccore", id= x, rettype="fasta_cds_aa", retmode="text")
    records = SeqIO.parse(handle, "fasta")
    for record in records:
        if record.description.__contains__('[pseudo=true]'):
            pass
        else:
            print(record.id)
            print(record.description)
            print(record.seq)
    #print(record.annotations["structured_comment"])
    handle.close()
    return record
"""

def plasmid_id():
    plasmids = pd.read_csv(unique_plasmids, index_col = 0, header = 0)
    plasmids = plasmids.rename(columns={"0": "Plasmid_ID"})
    id_list = plasmids['Plasmid_ID'].to_list()
    return id_list

def gb_retriever():
    id_list = plasmid_id()
    all_proteins = pd.DataFrame({'ID': [], 'P_Name': [], 'Sequence': []})
    for pl_id in id_list[:2]:
        all_info = search_Entrez(pl_id)
        all_proteins=all_proteins.append(all_info)
    return all_proteins

def df_to_fasta():
    df = gb_retriever()
    print(df)
    plsdb_orfeome = f'{tables}/plsdb_orfeome.fasta'
    ids = df['ID'].to_list()
    df['P_Name'] = df['P_Name'].apply(lambda x: x.replace(" ", "_"))
    print(df['P_Name'])
    prot_name = df['P_Name'].to_list()
    sequence = df['Sequence'].to_list()
    zipped = list(zip(ids, prot_name, sequence))
    records = []
    for record in zipped:
        id = record[0]
        description = record[1]
        sequence = record[2]
        the_record = SeqRecord(Seq(sequence), id = id, description = description)
        records.append(the_record)
        print(records)
    if not os.path.isfile(plsdb_orfeome) or os.stat(plsdb_orfeome).st_size == 0:
        SeqIO.write(records, plsdb_orfeome, "fasta")

df_to_fasta()
