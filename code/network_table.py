# -*- coding: utf-8 -*-
"""
Created on 09/12/2022 17:01

Author: Lucy
"""
### importing modules
import pandas as pd
import numpy as np
import os, re
from pathlib import Path

pd.set_option('display.max_columns', None)
### paths
# uncomment relevant path to OS
# Windows
path = r"C:\Users\Lucy\iCloudDrive\Documents\bengurion\Project students\Sivan_project"
# macOS
#path = r"/Users/lucyandrosiuk/Documents/bengurion/Project students/Sivan_project"
# Cluster
# path = r"/gpfs0/tals/projects/Analysis/Lucy_plasmidome/Plasmidome/CRISPR"

# working directories
visuals = f"{path}/visualisations"
tables = f"{path}/data_calculations"
resource = r"../res"
Path(visuals).mkdir(parents=True, exist_ok=True)

# working files
blast_results = f"{resource}/BLASTp_Database.zip"

def df_reader():
    df = pd.read_csv(blast_results)
    df = df[['qseqid', 'sseqid', 'ratio',  'plasmid family','spacer host taxonomy']]
    df = df.loc[df['qseqid'] != df['sseqid']]
    df[['Kingdom','Phylum','Class','Order','Family', 'Genus','else1', 'else2']] = df['spacer host taxonomy'].apply(lambda x: pd.Series(x.split(', ')))
    df = df.replace("'", '', regex=True)
    df = df.apply(lambda x: x.replace('[', '') if (isinstance(x, str)) else x)
    df = df.apply(lambda x: x.replace(']', '') if (isinstance(x, str)) else x)
    df.Family.fillna(df['plasmid family'], inplace = True)
    df = df[['qseqid', 'ratio',  'Family']]
    print(df)
    return df

def host():
    df = df_reader()
    plasmids = df['qseqid'].unique()
    print('There are %d unique plasmids' % df['qseqid'].nunique())
    families = df['Family'].unique()
    print('There are %d unique families' % df['Family'].nunique())


host()
