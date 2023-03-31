"""
Created on 30/03/2023

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

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

#pd.set_option('display.max_colwidth', None)
### paths
email = "androsiu@post.bgu.ac.il"  # Tell NCBI who you are
# uncomment relevant path to OS
# Windows
#path = r"C:\Users\Lucy\iCloudDrive\Documents\bengurion\Project students\Sivan_project"
# macOS
path = r"/Users/lucyandrosiuk/Documents/bengurion/Project students/Sivan_project"
# Cluster
#path = r"/gpfs0/tals/projects/Analysis/Lucy_plasmidome/Plasmidome/CRISPR"

# working directories
visuals = f"{path}/visualisations"
tables = f"{path}/data_calculations"
resource = f"{path}/res"
Path(visuals).mkdir(parents=True, exist_ok=True)
proteins_csv = f'{tables}/all_proteins_90.csv'
pseudo_proteins_csv = f'{tables}/pseudo_proteins_90.csv'

def work_files(cutoff):
    host_grades = f"{path}/cutoffs/id_{cutoff}/plasmids_grades_results_{cutoff}.csv"
    blank = f"{path}/cutoffs/id_{cutoff}/blank_results_{cutoff}.csv"
    match_update = f"{path}/match_update_{cutoff}.csv"
    all_info = f"{tables}/all_info_{cutoff}.csv"
    return host_grades, blank, match_update, all_info



def spacers_locations(cutoff):
    df_spacer = pd.read_csv(work_files(cutoff)[2], sep = ',', header = 0, index_col = 0)
    df_spacer = df_spacer[['qseqid', 'sseqid', 'length', 'ratio', 'qstart', 'qend', 'match']]
    df_range = pd.read_csv(work_files(cutoff)[3], sep = ',', header = 0, index_col = 0)
    df_range = df_range[['qseqid', 'level of difference', 'MOB']]
    # df_6 = df_range.loc[df_range['level of difference']==6]
    # plasmids = df_6['qseqid'].unique().tolist()
    print(df_range['qseqid'].unique().tolist())
    df_proteins = pd.read_csv(proteins_csv, sep = ',', header = 0, index_col = 0)
    #print(df_proteins.columns)

    print(df_proteins['qseqid'].nunique())
    df_proteins = df_proteins.rename(columns={"ID": "P_ID", 'Sequence': 'P_Sequence'})
    df_proteins = df_proteins[['P_ID', 'P_Name', 'P_start', 'P_end', 'qseqid']]

    df_pseudo_proteins = pd.read_csv(pseudo_proteins_csv, sep = ',', header = 0, index_col = 0)
    #print(df_pseudo_proteins.columns)

    print(df_pseudo_proteins['qseqid'].nunique())
    df_pseudo_proteins = df_pseudo_proteins.rename(columns={"ID": "Ps_ID", 'P_Name': 'Ps_Name', 'P_start':'Ps_start', 'P_end':'Ps_end', 'Sequence': 'Ps_Sequence'})
    df_pseudo_proteins = df_pseudo_proteins[['Ps_ID', 'Ps_Name', 'Ps_start', 'Ps_end','qseqid']]

    # create an empty list to store the matches
    all_matches = []
    ids = df_proteins['qseqid'].unique().tolist()

    # loop over each genome_id
    for genome_id in ids[150:160]:
        print(genome_id)
        # get the proteins, pseudo_proteins, and matches for the current genome_id
        proteins_df = df_proteins[df_proteins['qseqid'] == genome_id]
        #print(proteins_df)
        pseudo_proteins_df = df_pseudo_proteins[df_pseudo_proteins['qseqid'] == genome_id]
        #print(pseudo_proteins_df)
        matches_df = df_spacer[df_spacer['qseqid'] == genome_id]
        #print(matches_df)

        # loop over each protein row
        for protein_index, protein_row in proteins_df.iterrows():
            protein_start = protein_row['P_start']
            protein_end = protein_row['P_end']

            # loop over each pseudo_protein row
            for pseudo_protein_index, pseudo_protein_row in pseudo_proteins_df.iterrows():
                pseudo_protein_start = pseudo_protein_row['Ps_start']
                pseudo_protein_end = pseudo_protein_row['Ps_end']

                # loop over each match row
                for match_index, match_row in matches_df.iterrows():
                    match_start = match_row['qstart']
                    match_end = match_row['qend']

                    # check if the match falls within the protein coordinates
                    if (match_start >= protein_start) and (match_end <= protein_end):
                        # add the match to the list of matches that fall within proteins
                        all_matches.append({
                            'genome_id': genome_id,
                            'protein_id': protein_row['P_ID'],
                            'protein_name': protein_row['P_Name'],
                            'spacer_id': match_row['sseqid'],
                            'match_start': match_row['qstart'],
                            'match_end': match_row['qend'],
                            'within_protein': True,
                            'within_pseudo_protein': False,
                            'just_before_protein': False,
                            'just_after_protein': False
                        })
                    # check if the match falls within the pseudo_protein coordinates
                    elif (match_start >= pseudo_protein_start) and (match_end <= pseudo_protein_end):
                        # add the match to the list of matches that fall within pseudo_proteins
                        all_matches.append({
                            'genome_id': genome_id,
                            'protein_id': None,
                            'pseudo_protein_id': pseudo_protein_row['Ps_ID'],
                            'protein_name': pseudo_protein_row['Ps_Name'],
                            'spacer_id': match_row['sseqid'],
                            'match_start': match_row['qstart'],
                            'match_end': match_row['qend'],
                            'within_protein': False,
                            'within_pseudo_protein': True,
                            'just_before_protein': False,
                            'just_after_protein': False
                        })
                    # check if the match falls just before the protein
                    elif (match_end == protein_start):
                        all_matches.append({
                            'genome_id': genome_id,
                            'protein_id': protein_row['P_ID'],
                            'protein_name': protein_row['P_Name'],
                            'pseudo_protein_id': None,
                            'spacer_id': match_row['sseqid'],
                            'match_start': match_row['qstart'],
                            'match_end': match_row['qend'],
                            'within_protein': False,
                            'within_pseudo_protein': False,
                            'just_before_protein': True,
                            'just_after_protein': False
                        })
                    # check if the match falls just after the protein
                    elif (match_start == protein_end):
                        all_matches.append({
                            'genome_id': genome_id,
                            'protein_id': protein_row['P_ID'],
                            'protein_name': protein_row['P_Name'],
                            'pseudo_protein_id': None,
                            'spacer_id': match_row['sseqid'],
                            'match_start': match_row['qstart'],
                            'match_end': match_row['qend'],
                            'within_protein': False,
                            'within_pseudo_protein': False,
                            'just_before_protein': False,
                            'just_after_protein': True
                        })
                    else:
                        # add the match to the list of all matches, even if it doesn't fall within proteins or pseudo_proteins
                        all_matches.append({
                            'genome_id': genome_id,
                            'protein_id': None,
                            'protein_name': None,
                            'pseudo_protein_id': None,
                            'spacer_id': match_row['sseqid'],
                            'match_start': match_row['qstart'],
                            'match_end': match_row['qend'],
                            'within_protein': False,
                            'within_pseudo_protein': False,
                            'just_before_protein': False,
                            'just_after_protein': False
                        })
        print(all_matches)

    # create a new dataframe with all matches
    df_all_matches = pd.DataFrame(all_matches)
    print(df_all_matches.shape)
    print(df_all_matches)
    # write the results to a new CSV file
    df_all_matches.to_csv(f'{tables}/all_matches_{cutoff}.csv', index = False)

def cutoff():
    cutoff = [90, 95, 100]
    for i in cutoff:
        spacers_locations(i)

cutoff()