"""
Created on 26/03/2023

Author: Lucy Androsiuk
"""

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
resource = f"{path}/res"
Path(visuals).mkdir(parents=True, exist_ok=True)

# working files
copla_res = f"{tables}/results_tab.tsv"

def work_files(cutoff):
    host_grades = f"{path}/cutoffs/id_{cutoff}/plasmids_grades_results_{cutoff}.csv"
    blank = f"{path}/cutoffs/id_{cutoff}/blank_results_{cutoff}.csv"
    match_update = f"{path}/match_update_{cutoff}.csv"
    all_info = f"{tables}/all_info_{cutoff}.csv"
    return host_grades,blank,match_update, all_info


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
                p_start = int(p_loc[0])
                p_end = int(p_loc[-1])
                p_starts.append(p_start)
                p_ends.append(p_end)
                p = parts[-4]
                p_name = p[8:-2]
                p_names.append(p_name)
                sequence = str(record.seq)
                seqs.append(sequence)
        df = pd.DataFrame({'ID': ids, 'P_Name': p_names, 'P_start':p_starts, 'P_end':p_ends, 'Sequence': seqs})
        #print(df)
        print("############## Proteins for %s obtained ##############", x)
        #print(record.annotations["structured_comment"])
        handle.close()
        return df
    except Exception as e:
        print(e)


def plasmids(cutoff):
    df_spacer = pd.read_csv(work_files(cutoff)[2], sep = ',', header = 0, index_col = 0)
    df_spacer = df_spacer[['qseqid', 'sseqid', 'length', 'ratio', 'qstart', 'qend', 'match']]

    df_range = pd.read_csv(work_files(cutoff)[3], sep = ',', header = 0, index_col = 0)
    df_range = df_range[['qseqid', 'level of difference', 'MOB']]
    df_6 = df_range.loc[df_range['level of difference']==6]
    plasmids = df_6['qseqid'].unique().tolist()
    print(len(plasmids))
    df_proteins = pd.DataFrame()
    for plasmid in plasmids:
        df_proteins = df_proteins.append(search_Entrez(plasmid))
    df_proteins['qseqid'] = df_proteins['ID'].apply(lambda x: x.split('_', 1)[0])
    df_spacer = df_spacer.loc[df_spacer['qseqid'].isin(plasmids)]

    # merge the two dataframes based on the chromosome column
    merged_df = pd.merge(df_proteins, df_spacer, on = 'qseqid')
    # check if matches fall within the protein coordinates
    merged_df['match_within_protein'] = (merged_df['qstart'] >= merged_df['P_start']) & (
                merged_df['qend'] <= merged_df['P_end'])

    # filter the merged dataframe to show only matches that fall within proteins
    matches_within_proteins = merged_df.loc[merged_df['match_within_protein'] == True, :]
    print(matches_within_proteins.columns)
    matches_within_proteins_csv = f'{tables}/matches_within_proteins_{cutoff}.csv'
    if not os.path.isfile(matches_within_proteins_csv) or os.stat(matches_within_proteins_csv).st_size == 0:
        matches_within_proteins.to_csv(matches_within_proteins_csv, index = True)

def cutoff():
    cutoff = [90, 95, 100]
    for i in cutoff:
        plasmids(i)

cutoff()
