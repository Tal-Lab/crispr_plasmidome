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
import logging

#pd.set_option('display.max_colwidth', None)
### paths
email = "androsiu@post.bgu.ac.il"  # Tell NCBI who you are
# uncomment relevant path to OS
# Windows
#path = r"C:\Users\Lucy\iCloudDrive\Documents\bengurion\Project students\Sivan_project"
# macOS
#path = r"/Users/lucyandrosiuk/Documents/bengurion/Project students/Sivan_project"
# Cluster
path = r"/gpfs0/tals/projects/Analysis/Lucy_plasmidome/Plasmidome/CRISPR"

# working directories
visuals = f"{path}/visualisations"
tables = f"{path}/data_calculations"
Path(visuals).mkdir(parents=True, exist_ok=True)
Path(tables).mkdir(parents=True, exist_ok=True)

# working files
unique_plasmids = f"{path}/unique_plasmids.csv"
#unique_plasmids = r"../res/unique_plasmids.csv"

logging.basicConfig(filename=f'{path}/log.log', format = '%(asctime)s - %(levelname)s - %(message)s', level=logging.DEBUG)
logging.debug("I'm inside the python file")

def search_Entrez (x):
    '''getting orfeome for each spacer-barer id'''
    print("Entering search_Entrez function")
    Entrez.email = email  # Tell NCBI who you are
    Entrez.sleep_between_tries = 20  #Tell NCBI delay, in seconds, before retrying a request
    try:
        print("############## Obtaining proteins for %s ##############", x)
        handle = Entrez.efetch(db="nuccore", id= x, rettype="fasta_cds_aa", retmode="text")
        records = SeqIO.parse(handle, "fasta")
        ids = []
        p_names = []
        seqs = []
        p_starts = []
        p_ends = []
        for record in records:
            if record.description.__contains__('[pseudo=true]'):
                pass
            else:
                #print('Printing record')
                #print(record)
                id = record.id[4:]
                ids.append(id)
                description = record.description
                parts = description.split('[')
                p_loc_info = parts[-2]
                p_loc = re.findall(r'\d+', p_loc_info)
                p_start = p_loc[0]
                p_end = p_loc[-1]
                p_starts.append(p_start)
                p_ends.append(p_end)
                p = parts[-4]
                p_name = p[8:-2]
                p_names.append(p_name)
                sequence = str(record.seq)
                seqs.append(sequence)
        df = pd.DataFrame({'ID': ids, 'P_Name': p_names, 'P_start':p_starts, 'P_end':p_ends, 'Sequence': seqs})
        print(df)
        print("############## Proteins for %s obtained ##############", x)
        #print(record.annotations["structured_comment"])
        handle.close()
        return df
    except Exception as e:
        print(e)

def plasmid_id():
    logging.debug("Entering plasmid_id function")
    plasmids = pd.read_csv(unique_plasmids, index_col = 0, header = 0)
    plasmids = plasmids.rename(columns={"0": "Plasmid_ID"})
    id_list = plasmids['Plasmid_ID'].to_list()
    return id_list

def gb_retriever():
    logging.debug("Entering gb_retriever function")
    id_list = plasmid_id()
    all_proteins = pd.DataFrame({'ID': [], 'P_Name': [], 'Sequence': []})
    try:
        for pl_id in id_list:
            logging.info("### Working with plasmid %s ###", str(id_list.index(pl_id)))
            all_info = search_Entrez(pl_id)
            all_proteins=all_proteins.append(all_info)
    except Exception as e:
        print(e)
    return all_proteins

def df_to_fasta():
    logging.info("START")
    logging.debug("Entering df_to_fasta function")
    df = gb_retriever()
    plsdb_orfeome = f'{tables}/plsdb_orfeome.fasta'
    ids = df['ID'].to_list()
    df['P_Name'] = df['P_Name'].apply(lambda x: x.replace(" ", "_"))
    prot_name = df['P_Name'].to_list()
    sequence = df['Sequence'].to_list()
    zipped = list(zip(ids, prot_name, sequence))
    records = []
    try:
        logging.debug("############ CREATING FASTA WITH PLSDB ORFEOME #############")
        for record in zipped:
            id = record[0]
            description = record[1]
            sequence = record[2]
            the_record = SeqRecord(Seq(sequence), id = id, description = description)
            records.append(the_record)
        if not os.path.isfile(plsdb_orfeome) or os.stat(plsdb_orfeome).st_size == 0:
            SeqIO.write(records, plsdb_orfeome, "fasta")
            logging.debug("############ FASTA WITH PLSDB ORFEOME CREATED #############")
    except Exception as e:
        print(e)

df_to_fasta()
