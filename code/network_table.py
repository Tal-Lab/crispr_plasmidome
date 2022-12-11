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
import ast, random

pd.set_option('display.max_columns', None)

"""
1. https://stackoverflow.com/questions/110362/how-can-i-find-the-current-os-in-python
2. https://dev.to/jakewitcher/using-env-files-for-environment-variables-in-python-applications-55a1 depending on result of 1
"""

### paths
# uncomment relevant path to OS
# Windows
# path = r"C:\Users\Lucy\iCloudDrive\Documents\bengurion\Project students\Sivan_project"
# macOS
path = r"/Users/lucyandrosiuk/Documents/bengurion/Project students/Sivan_project"
# Cluster
# path = r"/gpfs0/tals/projects/Analysis/Lucy_plasmidome/Plasmidome/CRISPR"

# working directories
visuals = f"{path}/visualisations"
tables = f"{path}/data_calculations"
resource = r"../res"
Path(visuals).mkdir(parents=True, exist_ok=True)

# working files
blast_results = f"{resource}/BLASTp_Database.zip"

def colors_gen(x):
    print('Generating %d colors' % x)
    all_colors = []
    i = 0
    while i < x:
        r = random.randint(1, 255)
        g = random.randint(1, 255)
        b = random.randint(1, 255)
        rgb = f'{r}, {g}, {b}'
        if rgb != '0,152,152' and rgb not in all_colors:
            all_colors.append(rgb)
            i += 1
    return all_colors

def df_reader():
    df = pd.read_csv(blast_results)
    df = df[['qseqid', 'sseqid', 'ratio',  'plasmid family','spacer host taxonomy']]
    df = df.loc[df['qseqid'] != df['sseqid']]
    df['spacer host taxonomy'] = df['spacer host taxonomy'].apply(lambda x: ast.literal_eval(x))
    df[['Kingdom','Phylum','Class','Order','Family', 'Genus','else1', 'else2']] = pd.DataFrame(df['spacer host taxonomy'].tolist())
    df.Family.fillna(df['plasmid family'], inplace = True)
    df = df[['qseqid', 'Family']]
    df = df.drop_duplicates()
    print(df)
    network_csv = f'{tables}/plasmid_host_network3.csv'
    if not os.path.isfile(network_csv) or os.stat(network_csv).st_size == 0:
        df.to_csv(network_csv, index = False)
    pl = pd.DataFrame(df['qseqid'].dropna().unique(), columns=['id'])
    pl['type'] = 'plasmid'
    pl['colors'] = '0,152,152'
    pl['size'] = 20.0
    fam = pd.DataFrame(df['Family'].dropna().unique(), columns=['id'])
    fam['type'] = fam['id']
    fam['colors'] = colors_gen(len(fam))
    fam['size'] = 150.0
    df_type = pd.concat([pl,fam])
    print(df_type)
    nodes_csv = f'{tables}/nodes3.csv'
    if not os.path.isfile(nodes_csv) or os.stat(nodes_csv).st_size == 0:
        df_type.to_csv(nodes_csv, index = False)
    return df

df_reader()

