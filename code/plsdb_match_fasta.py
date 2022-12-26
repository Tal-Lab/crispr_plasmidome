# -*- coding: utf-8 -*-
"""
Created on 17/11/2022 14:33

Author: Lucy
"""
'''
Retrieving DNA sequence for unique plasmids with hits to spacers
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

### paths
# Cluster
path = r"/gpfs0/tals/projects/Analysis/Lucy_plasmidome/Plasmidome/CRISPR"

# working directories
visuals = f"{path}/visualisations"
tables = f"{path}/data_calculations"
resource = f"{path}/res"
Path(visuals).mkdir(parents=True, exist_ok=True)
Path(tables).mkdir(parents=True, exist_ok=True)

# working files
unique_plasmids = f"{resource}/unique_plasmids.csv"
plsdb_all =  r"/gpfs0/tals/projects/Analysis/sivan_project/plasmid_db/plsdb.fna"

logging.basicConfig(filename=f'{path}/log_uniquePlasmids.log', format = '%(asctime)s - %(levelname)s - %(message)s', level=logging.DEBUG)
logging.debug("I'm inside the python file")

def plasmid_id():
    logging.debug("Entering plasmid_id function")
    plasmids = pd.read_csv(unique_plasmids, index_col = 0, header = 0)
    plasmids = plasmids.rename(columns={"0": "Plasmid_ID"})
    id_list = plasmids['Plasmid_ID'].to_list()
    return id_list

def dna_retriever():
    unique = plasmid_id()
    names = []
    descs = []
    seqs = []
    with open(plsdb_all, 'r') as plsdb:
        for rec in SeqIO.parse(plsdb, 'fasta'):
            if rec.id in unique:
                name = rec.id
                desc = rec.description
                seq = str(rec.seq)
                names.append(name)
                descs.append(desc)
                seqs.append(seq)
            else:
                pass
    df = pd.DataFrame({'ID': names,'Description': descs, 'Sequence': seqs})
    return df

def df_to_fasta():
    logging.info("START")
    logging.debug("Entering df_to_fasta function")
    df = dna_retriever()
    plsdb_unique = f'{resource}/plsdb_unique.fasta'
    ids = df['ID'].to_list()
    descs = df['Description'].to_list()
    sequence = df['Sequence'].to_list()
    zipped = list(zip(ids, descs, sequence))
    records = []
    try:
        logging.debug("############ CREATING FASTA WITH PLSDB DNA #############")
        for record in zipped:
            id = record[0]
            description = record[1]
            sequence = record[2]
            the_record = SeqRecord(Seq(sequence), id = id, description = description)
            records.append(the_record)
        if not os.path.isfile(plsdb_unique) or os.stat(plsdb_unique).st_size == 0:
            SeqIO.write(records, plsdb_unique, "fasta")
            logging.debug("############ FASTA WITH PLSDB DNA CREATED #############")
    except Exception as e:
        logging.error(e)

df_to_fasta()


