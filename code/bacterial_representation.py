# -*- coding: utf-8 -*-
"""
Created on 19/07/2022 15:44

Author: Lucy
"""
'''
I need to write a function for hypergeomatric distribution here to detect whether families/classes are overrepresented in db.
'''

### importing modules
import pandas as pd
import numpy as np
from Bio import Entrez
from Bio import SeqIO
import os, re
from pathlib import Path

### paths
email = "androsiu@post.bgu.ac.il"  # Tell NCBI who you are
# uncomment relevant path to OS
# Windows
path = r"C:\Users\Lucy\iCloudDrive\Documents\bengurion\Project students\Sivan_project"
# macOS
#path = r"/Users/lucyandrosiuk/Documents/bengurion/Project students/Sivan_project"

# working directories
visuals = f"{path}/visualisations"
tables = f"{path}/data_calculations"
Path(visuals).mkdir(parents=True, exist_ok=True)

# working files
blast_results = r"../res/BLASTp_resultsRatio.zip"
plsdb_meta = r"../res/plsdb.zip"
spacers_meta = r"../res/spacer_seqName.fsa"

# managing metadata for spacers db
def search_Entrez (x):
    '''getting taxonomy information for each spacer-barer id'''
    Entrez.email = email  # Tell NCBI who you are
    Entrez.sleep_between_tries = 20 # Tell NCBI delay, in seconds, before retrying a request
    try:
        print("Entered search_Entrez function with %s" % x)
        handle = Entrez.efetch(db="nucleotide", id= x, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        print(record.annotations["taxonomy"], record.annotations["organism"])
        handle.close()
        return record.annotations["organism"], record.annotations["taxonomy"]
    except ValueError:
        return "Error"

def DF_from_fasta():
    '''making dataframe form spacers fasta with the bacterial ids to fetch taxonomy'''
    df_spacers = pd.DataFrame(columns = ["ID", "Sequence"])
    with open(spacers_meta, 'r') as spacers:
        for rec in SeqIO.parse(spacers, 'fasta'):
            name = rec.id
            seq = str(rec.seq)
            df=pd.DataFrame({'ID':[name],'Sequence':[seq]})
            df_spacers = df_spacers.append(df)
            #print(df_spacers)
            #df_spacers = df_spacers.assign(ID = df_spacers['ID'].str.split('+')).explode('ID')
            #print(df_spacers["ID"])
            #print("Getting taxonomy for %s" % df_spacers["ID"])
            #df_spacers["taxonomy"] = df_spacers["ID"].apply(search_Entrez)
            #print(df_spacers)
    df_spacers=df_spacers.assign(ID = df_spacers['ID'].str.split('+')).explode('ID')
    unique_id = df_spacers['ID'].unique()
    print(unique_id)
    taxon = []
    species = []
    taxon_df = pd.DataFrame(columns = ['ID','Taxonomy','Species'])
    taxon_csv = f'{tables}/taxon_efetch.csv'
    with open (taxon_csv, 'w') as csv_file:
        csv_file.write("ID,Taxonomy,Species\n")
        print("There are %d records" % len(unique_id))
        #ind=np.where(unique_id=="CP002775.1")
        #print("The index is %d" % ind)
        #ind = ind[0][0]
        #print("The index is %d" % ind)
        for id in unique_id:
            print(id)
            ind = np.where(unique_id == id)
            print("The index is %d" % ind)
            sp, tax = search_Entrez(id)
            df = pd.DataFrame({'ID':id, 'Taxonomy':[tax], 'Species':sp})
            df.to_csv(csv_file, index = False, header=False, mode = 'a')
            taxon.append((tax))
            species.append(sp)
    print(len(taxon))
    print(len(unique_id))
    print(len(species))
    taxon_df = pd.DataFrame({'ID':unique_id, 'Taxonomy':taxon, 'Species':species})
    print(taxon_df)
    df_final = pd.merge(df_spacers, taxon_df, on = 'ID', how = 'right')
    print(df_spacers.shape)
    final_csv = f'{tables}/blast_result.csv'
    df_final.to_csv(final_csv, index = False)
    print(df_final)
    print(df_final.shape)
    return df_final
DF_from_fasta()

# metadata for plasmid db
features_train = pd.read_csv(plsdb_meta, delimiter='\t', header=None,  error_bad_lines=False ,usecols=[1,28 ,30,32 , 34, 36, 38, 40])
#print(features_train)

# our results
plasmid_blast = pd.read_csv(blast_results, sep= ",",  header = 0)
#print(plasmid_blast)

# hypergeom distribution for each Genus, Family, Class, Phylum