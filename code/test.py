# -*- coding: utf-8 -*-
"""
Created on 02/04/2023 15:04

Author: Lucy
"""

import pandas as pd
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os, re


email = "androsiu@post.bgu.ac.il"  # Tell NCBI who you are

plasmids = ['AP019704.1', 'NZ_LT719075.1', 'NZ_CP050379.1']

def search_Entrez (x):
    '''getting orfeome for each spacer-barer id'''
    print("Entering search_Entrez function with plasmid ", x)
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
                print("!!!!!!!!!!!!!! PSEUDO")
                #print(record)
                id = record.id[4:]
                ids.append(id)
                description = record.description
                parts = description.split('[')
                #print(parts)
                p_loc_info = parts[-2]
                p_loc = re.findall(r'\d+', p_loc_info)
                p_start = int(p_loc[0])
                p_end = int(p_loc[-1])
                p_starts.append(p_start)
                p_ends.append(p_end)
                p_name = [match[8:-2] for match in parts if 'protein' in match]
                if p_name == []:
                    p_name = [match[11:-2] for match in parts if 'pseudogene' in match]
                print(p_name)
                p_names.append(p_name)
                sequence = str(record.seq)
                seqs.append(sequence)
            else:
                pass

            ''' 
            
            description = record.description
            parts = description.split('[')
            p_loc_info = parts[-2]
            p_loc = re.findall(r'\d+', p_loc_info)
            p_start = int(p_loc[0])
            p_end = int(p_loc[-1])
            p_starts.append(p_start)
            p_ends.append(p_end)
            p = parts[-4]
            p_name = p[8:-2]
            p_names.append(p_name)
            sequence = str(record.seq)
            seqs.append(sequence)
            '''
        df = pd.DataFrame({'ID': ids, 'P_Name': p_names, 'P_start':p_starts, 'P_end':p_ends, 'Sequence': seqs})
        #print(df)
        print("############## Proteins for %s obtained ##############", x)
        #print(record.annotations["structured_comment"])
        handle.close()
        return df
    except Exception as e:
        print(e)

for plasmid in plasmids:
    search_Entrez(plasmid)