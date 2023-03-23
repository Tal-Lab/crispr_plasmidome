# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 19:20:44 2022

@author: sivan
"""
 
import pandas as pd
import numpy as np
import sys
from multiprocessing import Pool
from time import process_time
import itertools
from functools import partial

#plasmid_blast = pd.read_csv (sys.argv[1], sep= ",",  header = 0)
plasmid_blast = pd.read_csv (r"C:\Users\Lucy\Dev\crispr_plasmids\res\BLASTp_DataBase3.zip", sep= ",",  header = 0)
#remove all rows where qseqid == sseqid
concat_blast = plasmid_blast[plasmid_blast.qseqid != plasmid_blast.sseqid]

def Convert(string):
    string_clean = string [1:-1].replace("'", "")    
    li = list(string_clean.split(","))
    return li

concat_blast['spacer host taxonomy'] = concat_blast['spacer host taxonomy'].apply(Convert)

#remove groups that are not relevent to the taxonomy hierarchy 
def remove_groups(taxonomy):
    try:
        if 'group' in taxonomy[5]:
            taxonomy.pop(5)
    except IndexError:
        return taxonomy
    return taxonomy

concat_blast['spacer host taxonomy'] = concat_blast['spacer host taxonomy'].apply(remove_groups)

#setting ratio, identity and mismatch cutoff
ratio_cutoff = 0.8 #constant
mismatch = 2 #constant

def pident_cutoff(identity_cutoff, df):
    df = df.loc[df['ratio'] >= ratio_cutoff]
    df = df.loc[df['mismatch'] <= mismatch]
    df = df.loc[df['pident'] >= identity_cutoff]
    print('printing after %s cutoff' % identity_cutoff)
    print(df.shape)
    plasmids_name = f'plasmids_{identity_cutoff}.csv'
    #plasmids_name = f"C:/Users/Lucy/iCloudDrive/Documents/bengurion/Project students/Sivan_project/match_update_{identity_cutoff}.csv"
    print(plasmids_name)
    df.to_csv(plasmids_name)
    #making a list of qseqid
    list_of_qseqid = df['qseqid'].unique()
    dict = {'Plasmid':list_of_qseqid}
    df_list = pd.DataFrame(dict)
    list_name = f'unique_plasmids_{identity_cutoff}.txt'
    df_list.to_csv((list_name))
    return df, list_of_qseqid

def host_range(cutoff, df):
    # creating new dataframe
    df = pident_cutoff(cutoff, df)[0]
    df_of_qseqid = pd.DataFrame(df['qseqid'].unique())
    df_of_qseqid['level of difference'] = ""
    df_of_qseqid.rename(columns = {0: 'qseqid'}, inplace = True)
    # go over qseqid list and create a dataframe from all result of the specific id in blast
    list_of_qseqid = pident_cutoff(cutoff, df)[1]
    for i in list_of_qseqid:
        t1_start = process_time()
        rslt_df = concat_blast.loc[concat_blast['qseqid'] == i]
        rslt_df = rslt_df.reset_index()
        rslt_df['level of difference'] = ""
        # go over rows of the temporary datafrme of the specific id and check if there are results for different hosts(taxonomy), if so check the level of difference
        for e in rslt_df.index:
           try:
                if len(rslt_df['spacer host taxonomy'][0]) != len(rslt_df['spacer host taxonomy'][e]):
                    if rslt_df['spacer host taxonomy'][0][1] != rslt_df['spacer host taxonomy'][e][1]:
                        rslt_df.loc[e, 'level of difference'] = "Phylum"
                    elif rslt_df['spacer host taxonomy'][0][2] != rslt_df['spacer host taxonomy'][e][2]:
                        rslt_df.loc[e, 'level of difference'] = "Class"
                    elif rslt_df['spacer host taxonomy'][0][3] != rslt_df['spacer host taxonomy'][e][3]:
                        rslt_df.loc[e, 'level of difference'] = 'Order'
                    elif rslt_df['spacer host taxonomy'][0][4] != rslt_df['spacer host taxonomy'][e][4]:
                        rslt_df.loc[e, 'level of difference'] = 'Family'
                    elif rslt_df['spacer host taxonomy'][0][5] != rslt_df['spacer host taxonomy'][e][5]:
                        rslt_df.loc[e, 'level of difference'] = 'Genus'
                    elif rslt_df['spacer host species'][0] != rslt_df['spacer host species'][e]:
                        rslt_df.loc[e, 'level of difference'] = 'Species'
                else:
                    if rslt_df['spacer host taxonomy'][0] != rslt_df['spacer host taxonomy'][e]:
                        if rslt_df['spacer host taxonomy'][0][1] != rslt_df['spacer host taxonomy'][e][1]:
                            rslt_df.loc[e, 'level of difference'] = "Phylum"
                        elif rslt_df['spacer host taxonomy'][0][2] != rslt_df['spacer host taxonomy'][e][2]:
                            rslt_df.loc[e, 'level of difference'] = "Class"
                        elif rslt_df['spacer host taxonomy'][0][3] != rslt_df['spacer host taxonomy'][e][3]:
                            rslt_df.loc[e, 'level of difference'] = 'Order'
                        elif rslt_df['spacer host taxonomy'][0][4] != rslt_df['spacer host taxonomy'][e][4]:
                            rslt_df.loc[e, 'level of difference'] = 'Family'
                        elif rslt_df['spacer host taxonomy'][0][5] != rslt_df['spacer host taxonomy'][e][5]:
                            rslt_df.loc[e, 'level of difference'] = 'Genus'
                        elif rslt_df['spacer host species'][0] != rslt_df['spacer host species'][e]:
                            rslt_df.loc[e, 'level of difference'] = 'Species'
                    else:
                        if rslt_df['spacer host species'][0] != rslt_df['spacer host species'][e]:
                            rslt_df.loc[e, 'level of difference'] = 'Species'

           except IndexError:
                if rslt_df['spacer host species'][0] != rslt_df['spacer host species'][e]:
                    rslt_df.loc[e, 'level of difference'] = 'Species'
        # check the most far taxonomy difference
        if 'Phylum' in rslt_df['level of difference'].unique():
            df_of_qseqid['level of difference'][df_of_qseqid['qseqid'] == i] = "6"
        elif 'Class' in rslt_df['level of difference'].unique():
            df_of_qseqid['level of difference'][df_of_qseqid['qseqid'] == i] = "5"
        elif "Order" in rslt_df['level of difference'].unique():
            df_of_qseqid['level of difference'][df_of_qseqid['qseqid'] == i] = "4"
        elif "Family" in rslt_df['level of difference'].unique():
            df_of_qseqid['level of difference'][df_of_qseqid['qseqid'] == i] = "3"
        elif "Genus" in rslt_df['level of difference'].unique():
            df_of_qseqid['level of difference'][df_of_qseqid['qseqid'] == i] = "2"
        elif "Species" in rslt_df['level of difference'].unique():
            df_of_qseqid['level of difference'][df_of_qseqid['qseqid'] == i] = "1"
    t1_stop = process_time()
    print("Elapsed time in seconds:", t1_stop - t1_start)
    print('printing with level of difference')
    print(df_of_qseqid.shape)
    grades_name = f'plasmids_grades_results_{cutoff}.csv'
    print(grades_name)
    df_of_qseqid.to_csv(grades_name)
    # create datafrme to check the resone for rows with no results
    # if their is one hit for the qseqid, the resukt will be 1
    # if all hits are for the same host, the result will be equal
    t2_start = process_time()
    tut = df_of_qseqid.loc[df_of_qseqid['level of difference'] == ""]
    tut = tut.reset_index()
    for d in tut['qseqid']:
        count_rows = concat_blast.loc[concat_blast['qseqid'] == d]
        count_rows = count_rows.reset_index()
        count_rows['equal'] = ""
        if len(count_rows.index) == 1:
            tut['level of difference'][tut['qseqid'] == d] = 1
        else:
            for s in count_rows.index:
                if len(count_rows['spacer host taxonomy'][0]) == len(count_rows['spacer host taxonomy'][s]):
                    if count_rows['spacer host taxonomy'][s] == count_rows['spacer host taxonomy'][0]:
                        count_rows['equal'][s] = 'yes'
                    else:
                        count_rows['equal'][s] = 'no'
                else:
                    count_rows['equal'][s] = 'no'
        if 'no' in count_rows['equal'].unique():
            tut['level of difference'][tut['qseqid'] == d] = 'check'
        elif 'yes' in count_rows['equal'].unique():
            tut['level of difference'][tut['qseqid'] == d] = 'equal'
    t2_stop = process_time()
    print("Elapsed time during the whole program in seconds:", t2_stop - t2_start)
    print('printing blanks')
    print(tut.shape)
    empty_name = f'blank_results_{cutoff}.csv'
    tut.to_csv(empty_name)

def parallel_runs(data_list):
    with Pool(4) as p:
        pident_cutoff_x = partial(host_range, df = concat_blast)
        result_list = p.map(pident_cutoff_x, data_list)
        print(result_list)
'''
if __name__ == '__main__':
    cutoff = [90,95,100]
    t_start = process_time()
    parallel_runs(cutoff)
    t_stop = process_time()
    print("Elapsed time during processing in seconds:", t_stop - t_start)
'''

def spacer_statistics(df):
    cutoff = [90, 95, 100]
    for i in cutoff:
        df_i = pident_cutoff(i, df)[0]
        print(df_i.columns)
        df_un = df_i.drop_duplicates(['qseqid','sseqid'])
        df_group = df_un.groupby('qseqid')['sseqid'].count()
        stat = f"C:/Users/Lucy/iCloudDrive/Documents/bengurion/Project students/Sivan_project/data_calculations/spacer_statisitcs_{i}.csv"
        df_group.to_csv(stat)

def species_statistics(df):
    cutoff = [90, 95, 100]
    for i in cutoff:
        df_i = pident_cutoff(i, df)[0]
        print(df_i.columns)
        df_un = df_i.drop_duplicates(['qseqid','spacer host species'])
        df_group = df_un.groupby('qseqid')['spacer host species'].count()
        stat = f"C:/Users/Lucy/iCloudDrive/Documents/bengurion/Project students/Sivan_project/data_calculations/species_statisitcs_{i}.csv"
        df_group.to_csv(stat)

def match_update():
    df =  pd.read_csv (r"C:\Users\Lucy\iCloudDrive\Documents\bengurion\Project students\Sivan_project\match_update.zip", sep= ",",  header = 0)
    print(df)
    cutoff = [90, 95, 100]
    for i in cutoff:
        pident_cutoff(i, df)

#match_update()



#spacer_statistics(concat_blast)
#species_statistics(concat_blast)

