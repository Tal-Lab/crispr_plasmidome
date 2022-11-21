# -*- coding: utf-8 -*-
"""
Created on 19/07/2022 15:44

Author: Lucy
"""
'''
Function for hypergeomatric distribution to detect whether families/classes are overrepresented in db.
'''

### importing modules
import pandas as pd
import numpy as np
from Bio import Entrez
from Bio import SeqIO
import os, re
from pathlib import Path
from scipy.stats import hypergeom

pd.set_option('display.max_columns', None)
### paths
email = "androsiu@post.bgu.ac.il"  # Tell NCBI who you are
# uncomment relevant path to OS
# Windows
path = r"C:\Users\Lucy\iCloudDrive\Documents\bengurion\Project students\Sivan_project"
# macOS
#path = r"/Users/lucyandrosiuk/Documents/bengurion/Project students/Sivan_project"
# Cluster
#path = r"/gpfs0/tals/projects/Analysis/Lucy_plasmidome/Plasmidome/CRISPR"

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
        handle.close()
        return record.annotations["organism"], record.annotations["taxonomy"]
    except ValueError:
        return "Error"

def DF_from_fasta():
    '''making dataframe form spacers fasta with the bacterial ids to fetch taxonomy'''
    df_spacers = pd.DataFrame(columns = ["ID", "Sequence"]) # creating empty dataframe for spacers' host-ids
    with open(spacers_meta, 'r') as spacers:
        for rec in SeqIO.parse(spacers, 'fasta'):
            name = rec.id
            seq = str(rec.seq)
            df=pd.DataFrame({'ID':[name],'Sequence':[seq]})
            df_spacers = df_spacers.append(df)
    df_spacers=df_spacers.assign(ID = df_spacers['ID'].str.split('+')).explode('ID')
    unique_id = df_spacers['ID'].unique() #obtaining list of unique host ids
    taxon = []
    species = []
    taxon_csv = f'{tables}/taxon_efetch.csv'
    with open (taxon_csv, 'w') as csv_file:
        csv_file.write("ID,Taxonomy,Species\n")
        print("There are %d records" % len(unique_id))
        for id in unique_id[:3]:
            ind = np.where(unique_id == id) # retrieving id's index
            print("The index is %d" % ind)
            sp, tax = search_Entrez(id) # retrieving species and taxonomy info for each id
            #df = pd.DataFrame({'ID':id, 'Taxonomy':[tax], 'Species':sp}) # creating dataframe with spacers' host lineage per id
            #df.to_csv(csv_file, index = False, header=False, mode = 'a')
            taxon.append((tax))
            species.append(sp)
    taxon_df = pd.DataFrame({'ID':unique_id[:3], 'Taxonomy':taxon, 'Species':species}) # creating dataframe with spacers' host lineage
    df_final = pd.merge(df_spacers, taxon_df, on = 'ID', how = 'right')
    df_final['Genus'] = df_final['Taxonomy'].apply(lambda x: x[5])
    df_final['Family'] = df_final['Taxonomy'].apply(lambda x: x[4])
    df_final['Order'] = df_final['Taxonomy'].apply(lambda x: x[3])
    df_final['Class'] = df_final['Taxonomy'].apply(lambda x: x[2])
    df_final['Phylum'] = df_final['Taxonomy'].apply(lambda x: x[1])
    df_final['Kingdom'] = df_final['Taxonomy'].apply(lambda x: x[0])
    df_final.drop('Taxonomy', axis=1, inplace = True)
    print(df_final)
    final_csv = f'{tables}/.csv'
    #df_final.to_csv(final_csv, index = False)
    return df_final
#DF_from_fasta()

def plsdb():
    # metadata for plasmid db
    plsdb_hosts = pd.read_csv(plsdb_meta, delimiter='\t', header=0,  error_bad_lines=False ,usecols=[1,28 ,30,32 , 34, 36, 38, 40])
    print(plsdb_hosts)
    return plsdb_hosts


def blast_res():
    # our results
    plasmid_blast = pd.read_csv(blast_results, sep= ",",  header = 0)
    df_spacers = DF_from_fasta()
    df_plasmids = plsdb()
    print(plasmid_blast)
    blast_df = pd.merge(plasmid_blast, df_plasmids, left_on = 'qseqid' , right_on =  'ACC_NUCCORE', how = 'left')
    blast_df = pd.merge(blast_df, df_spacers, left_on = 'stitle' , right_on = 'ID', how = 'left')
    print(blast_df)
    return blast_df
df_blast = blast_res()

# hypergeom distribution for each Genus, Family, Class, Phylum
def prob_func(x, level, popul, df_popul, df_blast):
    ''' Function to calculate COG categories statistics. \
    Function gets COG category (x), taxonomy level (level), number of plasmids/hosts, \
    with any host (hosts), spacer or plasmid database (df), blast_result (df_ps).'''
    if popul == 'plasmids':
        level_name = 'taxon_'+ level + '_name'
    elif popul == 'spacers':
        level_name = level.capitalize()
    else:
        print('something went wrong with level names')
    # getting N - number of spacers/plasmids (population size)
    N = len(df_popul)
    print(N)

    # getting A - number of plasmids/spacers with particular host x from spacer/plasmid (df_popul)
    df_popul = df_popul[df_popul[level_name] == x]
    A = len(df_popul)
    print(A)

    # getting n - number of hits
    n = len(df_blast)
    print(n)

    # getting k - number of plasmids/spacers with particular host x from blast_result databases
    df_blast = df_blast[df_blast[level_name] == x]
    k = len (df_blast)
    print(k)

    print(N, A, n, k)
    rv = hypergeom(N, A, n)
    pmf_level = rv.pmf(k)
    print(pmf_level)
    pval = hypergeom.sf(k-1, N, A, n)
    print(pval)

    level_freq = (A/N)*100
    level_blast_freq = (k/n)*100
    #print (pmf_cog)
    print('The probability of getting %d %s with host %s at %s level out of %d hits is %s.' % (k, popul, x, level, n, "{:.2e}".format(pmf_level)))
    print('The probability of getting %d or more %s with host %s at %s level out of %d hits is %s.' % (k, popul, x, level, n, "{:.2e}".format(pval)))
    print('The frequency of %s with host %s at %s level in %s-database: (%d/%d)*100=%f' % (popul, x, level, popul, A, N, round(level_freq,2)))
    print('The frequency of %s with host %s at %s in our blast results: (%d/%d)*100=%f' % (popul, x, level, k, n, round(level_blast_freq, 2)))
    #not working
    prob = hypergeom.pmf(N, A, n, k)
    print(prob)
    return pmf_level

def getting_levels(subj):
    if subj == 'plasmids':
        un_genus = df_blast['taxon_genus_name'].unique()
        un_family = df_blast['taxon_family_name'].unique()
        un_class = df_blast['taxon_class_name'].unique()
        un_phylum = df_blast['taxon_phylum_name'].unique()
    elif subj == 'spacers':
        un_genus = df_blast['Genus'].unique()
        un_family = df_blast['Family'].unique()
        un_class = df_blast['Class'].unique()
        un_phylum = df_blast['Phylum'].unique()
    else:
        print("something went wrong")
    return [un_genus, un_family, un_class, un_phylum]

def getting_probability(subj, df_popul):
    level_lists = getting_levels(subj)
    levels = ['genus', 'family', 'class', 'phylum']
    for level in levels:
        id = levels.index(level)
        level_list = level_lists[id]
        for taxa in level_list:
            prob_func(taxa, level, subj, df_popul, df_blast)

getting_probability('plasmids', plsdb())

