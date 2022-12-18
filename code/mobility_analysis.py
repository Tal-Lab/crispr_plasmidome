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
path = r"C:\Users\Lucy\iCloudDrive\Documents\bengurion\Project students\Sivan_project"
# macOS
#path = r"/Users/lucyandrosiuk/Documents/bengurion/Project students/Sivan_project"

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

def read_copla():
    df = pd.read_csv(copla_res, delimiter = '\t', header = None)
    df = df[[0, 1]]
    colnames = ['query', 'MOB']
    df.columns = colnames
    df['Plasmid'] = df['query'].apply(extract_name)
    df = df[['Plasmid', 'MOB']]
    df = df.drop_duplicates()
    return df

def mobile_copla():
    df = read_copla()
    df = df.drop_duplicates('Plasmid')
    plasmids = df['Plasmid'].tolist()
    return plasmids

def mob_grades():
    df = pd.read_csv(host_grades, index_col = 0)
    mob_plasmids = mobile_copla()
    df['MOB'] = df['qseqid'].apply(lambda x: 'MOB+' if x in mob_plasmids else 'MOB-')
    df['level of difference'] = df['level of difference'].fillna(1)
    #print(df)
    mobile = df.loc[df['MOB'] == 'MOB+']
    nonmobile = df.loc[df['MOB'] == 'MOB-']
    mob_csv = f'{tables}/mobile_grades.csv'
    if not os.path.isfile(mob_csv) or os.stat(mob_csv).st_size == 0:
        mobile.to_csv(mob_csv, index = True)
    nonmob_csv = f'{tables}/non_mobile_grades.csv'
    if not os.path.isfile(nonmob_csv) or os.stat(nonmob_csv).st_size == 0:
        nonmobile.to_csv(nonmob_csv, index = True)
    return df

def selectFromTable(df, criterias):
    for index, criteria in enumerate(criterias):
        ### checking for each criteria until the condition is True and returning best match
        filtered = df.loc[df['match'] == criteria]
        if len(filtered) > 0:
            return filtered

def table_all_info():
    ''' This function generates general table containing all information about
        Plasmid, plsdb host lineage, host grade, mobility'''
    ### reading table with matches
    df = pd.read_csv(r'/Users/lucyandrosiuk/Documents/bengurion/Project students/Sivan_project/match_update.csv',header=0, index_col = [0,1])
    df = df[['qseqid', 'plasmid species', 'plasmid genus', 'plasmid family', 'plasmid order', 'plasmid class','plasmid phylum', 'match']]
    df.drop_duplicates(inplace = True)

    ### group table by Plasmids to select only best match - the closest to PLSDB
    grouped = df.groupby('qseqid')
    # creating empty df for best matches
    best_match = pd.DataFrame(columns = ['qseqid', 'plasmid species', 'plasmid genus', 'plasmid family', 'plasmid order', 'plasmid class','plasmid phylum', 'match'])
    criteria_list = ['Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Superkingdom']
    for name, group in grouped:
        ### getting best match in each group (for each Plasmid) by the lowest match level
        matched = selectFromTable(group, criteria_list)
        # appending result to best match df
        best_match = best_match.append(matched)

    ### obtaining mob+ and mob- dataframes from mob_grades(), filling NAs with 1, and concatenating
    mob = mob_grades()

    ### merging plasmid matches df with mobility df
    table = best_match.merge(mob, how='left', on = 'qseqid')
    ### writing resulting dataframe into csv
    all_info = f'{tables}/all_info.csv'
    if not os.path.isfile(all_info) or os.stat(all_info).st_size == 0:
        table.to_csv(all_info, index = True)

#table_all_info()