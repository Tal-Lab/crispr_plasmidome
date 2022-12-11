"""
Created on 11/12/2022

Author: Lucy Androsiuk
"""

### importing modules
import pandas as pd
import numpy as np
import os, re
from pathlib import Path

# uncomment relevant path to OS
# Windows
#path = r"C:\Users\Lucy\iCloudDrive\Documents\bengurion\Project students\Sivan_project"
# macOS
path = r"/Users/lucyandrosiuk/Documents/bengurion/Project students/Sivan_project"

# working directories
visuals = f"{path}/visualisations"
tables = f"{path}/data_calculations"
resource = f"{path}/res"
Path(visuals).mkdir(parents=True, exist_ok=True)

# working files
copla_res = f"{tables}/results_tab.tsv"
host_grades = f"{path}/plasmids_grades_rsults_update2/plasmids_grades_rsults_update2.csv"

def extract_name(x):
    try:
        pl_name = re.search(r'.+(?=_prot)', x).group(0)
        return pl_name
    except:
        return x

def mobile_copla():
    df = pd.read_csv(copla_res, delimiter = '\t', header = None)
    df=df[[0,1]]
    colnames = ['query', 'MOB']
    df.columns = colnames
    df['Plasmid'] = df['query'].apply(extract_name)
    df = df[['Plasmid', 'MOB']]
    df = df.drop_duplicates()
    df = df.drop_duplicates('Plasmid')
    plasmids = df['Plasmid'].tolist()
    return plasmids

def mob_grades():
    df = pd.read_csv(host_grades, index_col = 0)
    mob_plasmids = mobile_copla()
    df['MOB'] = df['qseqid'].apply(lambda x: 'MOB+' if x in mob_plasmids else 'MOB-')
    print(df)
    mobile = df.loc[df['MOB'] == 'MOB+']
    nonmobile = df.loc[df['MOB'] == 'MOB-']
    mob_csv = f'{tables}/mobile_grades.csv'
    if not os.path.isfile(mob_csv) or os.stat(mob_csv).st_size == 0:
        mobile.to_csv(mob_csv, index = True)
    nonmob_csv = f'{tables}/non_mobile_grades.csv'
    if not os.path.isfile(nonmob_csv) or os.stat(nonmob_csv).st_size == 0:
        nonmobile.to_csv(nonmob_csv, index = True)

mob_grades()