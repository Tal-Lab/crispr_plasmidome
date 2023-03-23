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

def work_files(cutoff):
    host_grades = f"{tables}/plasmids_grades_results_{cutoff}.csv"
    blank = f"{tables}/blank_results_{cutoff}.csv"
    match_update = f"{path}/match_update_{cutoff}.csv"
    return host_grades,blank,match_update

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

def blank_grades(df, cutoff):
    df_blank = pd.read_csv(work_files(cutoff)[1], header = 0, index_col = 0)
    df_blank.drop(['index'], axis = 1, inplace = True)
    df_1 = df_blank.loc[df_blank['level of difference'] != 'equal']
    plasmids_to_remove = df_1['qseqid'].unique().tolist()
    #print(plasmids_to_remove)
    df['qseqid'] = df['qseqid'].apply(lambda x: x if x not in plasmids_to_remove else 'remove')
    df = df.loc[df['qseqid'] != 'remove']
    df['level of difference'].fillna(1, inplace = True)
    #df.Family.fillna(df['plasmid family'], inplace = True)
    return df

def mob_grades(cutoff):
    df = pd.read_csv(work_files(cutoff)[0], index_col = 0)
    no_blank = blank_grades(df,cutoff)
    print(no_blank)
    mob_plasmids = mobile_copla()
    #no_blank['MOB'] = no_blank['qseqid'].apply(lambda x: 'MOB+' if x in mob_plasmids else 'MOB-')
    no_blank.loc[no_blank['qseqid'].isin(mob_plasmids), 'MOB'] = 'MOB+'
    no_blank.loc[~no_blank['qseqid'].isin(mob_plasmids), 'MOB'] = 'MOB-'
    #no_blank['level of difference'] = no_blank['level of difference'].fillna(1)
    #print(df)
    mobile =  no_blank.loc[no_blank['MOB'] == 'MOB+']
    nonmobile = no_blank.loc[no_blank['MOB'] == 'MOB-']
    mob_csv = f'{tables}/mobile_grades_{cutoff}.csv'
    if not os.path.isfile(mob_csv) or os.stat(mob_csv).st_size == 0:
        mobile.to_csv(mob_csv, index = True)
    nonmob_csv = f'{tables}/non_mobile_grades_{cutoff}.csv'
    if not os.path.isfile(nonmob_csv) or os.stat(nonmob_csv).st_size == 0:
        nonmobile.to_csv(nonmob_csv, index = True)
    return no_blank

def selectFromTable(df, criterias):
    for index, criteria in enumerate(criterias):
        ### checking for each criteria until the condition is True and returning best match
        filtered = df.loc[df['match'] == criteria]
        if len(filtered) > 0:
            return filtered

def table_all_info(cutoff):
    ''' This function generates general table containing all information about
        Plasmid, plsdb host lineage, host grade, mobility'''
    ### reading table with matches
    df = pd.read_csv(work_files(cutoff)[2],header=0, index_col = 0)
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
    mob = mob_grades(cutoff)
    plasmids = mob['qseqid'].unique().tolist() # obtaining plasmids
    ### merging plasmid matches df with mobility df
    best_match = best_match.loc[best_match['qseqid'].isin(plasmids)]
    table = best_match.merge(mob, how='left', on = 'qseqid')
    ### writing resulting dataframe into csv
    all_info = f'{tables}/all_info_{cutoff}.csv'
    if not os.path.isfile(all_info) or os.stat(all_info).st_size == 0:
        table.to_csv(all_info, index = True)
'''
cutoff = [90, 95, 100]
for i in cutoff:
    table_all_info(i)
'''