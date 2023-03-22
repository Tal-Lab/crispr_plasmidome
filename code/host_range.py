# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 19:20:44 2022

@author: sivan
"""
 
import pandas as pd
import numpy as np
import sys

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
        print(taxonomy)
    return taxonomy

concat_blast['spacer host taxonomy'] = concat_blast['spacer host taxonomy'].apply(remove_groups)


#setting ratio cutoff
ratio_cutoff = 0.95
identity_cutoff = 95
print('printing before cutoff')
print(concat_blast.shape)
concat_blast = concat_blast.loc[concat_blast['ratio'] >= ratio_cutoff]
print('printing after cutoff')
print(concat_blast.shape)

concat_blast = concat_blast.loc[concat_blast['mismatch'] <= 2]
print('printing after cutoff')
print(concat_blast.shape)

concat_blast = concat_blast.loc[concat_blast['pident'] >= identity_cutoff]

print('printing after cutoff')
print(concat_blast.shape)
concat_blast = concat_blast.loc[concat_blast['pident'] >= 100]
print('printing after cutoff')
print(concat_blast.shape)

'''
#making a list of qseqid
list_of_qseqid = concat_blast ['qseqid'].unique()


#creating new dataframe
df_of_qseqid = pd.DataFrame(concat_blast['qseqid'].unique())
df_of_qseqid['level of difference']= ""
df_of_qseqid.rename(columns = {0:'qseqid'}, inplace = True)



#go over qseqid list and create a dataframe from all result of the specific id in blast
for i in list_of_qseqid:
    print('i am going over qseqid list and create a dataframe from all result of the specific id in blast')
    rslt_df = concat_blast.loc[concat_blast['qseqid'] == i]
    rslt_df = rslt_df.reset_index()
    rslt_df['level of difference'] = ""
    #go over rows of the temporary datafrme of the specific id and check if there are results for different hosts(taxonomy), if so check the level of difference
    for e in rslt_df.index:
        print('i am going over rows of the temporary datafrme of the specific id and check if there are results for different hosts(taxonomy), if so check the level of difference')
        try:
            if len(rslt_df['spacer host taxonomy'][0]) != len(rslt_df['spacer host taxonomy'][e]):
                if rslt_df['spacer host taxonomy'][0][1] != rslt_df['spacer host taxonomy'][e][1]:
                    rslt_df.loc[e,'level of difference'] = "Phylum"
                elif rslt_df['spacer host taxonomy'][0][2] != rslt_df['spacer host taxonomy'][e][2]:
                    rslt_df.loc[e,'level of difference'] = "Class"
                elif rslt_df['spacer host taxonomy'][0][3] != rslt_df['spacer host taxonomy'][e][3]:
                    rslt_df.loc[e,'level of difference'] = 'Order'
                elif rslt_df['spacer host taxonomy'][0][4] != rslt_df['spacer host taxonomy'][e][4]:
                    rslt_df.loc[e,'level of difference'] = 'Family'
                elif rslt_df['spacer host taxonomy'][0][5] != rslt_df['spacer host taxonomy'][e][5]:
                    rslt_df.loc[e,'level of difference'] = 'Genus'   
                elif rslt_df['spacer host species'][0] != rslt_df['spacer host species'][e]:
                    rslt_df.loc[e,'level of difference'] = 'Species'
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
                rslt_df.loc[e,'level of difference'] = 'Species' 
    #check the most far taxonomy difference
    print('i am checking the most far taxonomy difference')
    if 'Phylum' in rslt_df['level of difference'].unique():
        df_of_qseqid['level of difference'][df_of_qseqid['qseqid'] == i ]= "6"
    elif 'Class' in rslt_df['level of difference'].unique():
        df_of_qseqid['level of difference'][df_of_qseqid['qseqid'] == i ] = "5"
    elif "Order" in rslt_df['level of difference'].unique():
        df_of_qseqid['level of difference'] [df_of_qseqid['qseqid'] == i ] ="4"
    elif "Family" in rslt_df['level of difference'].unique():
        df_of_qseqid['level of difference'][df_of_qseqid['qseqid'] == i] = "3"
    elif "Genus" in rslt_df['level of difference'].unique():
        df_of_qseqid ['level of difference'][df_of_qseqid['qseqid'] == i]= "2"
    elif "Species" in rslt_df['level of difference'].unique():
        df_of_qseqid ['level of difference'][df_of_qseqid['qseqid'] == i] = "1"
   
        
#create datafrme to check the resone for rows with no results
#if their is one hit for the qseqid, the resukt will be 1
#if all hits are for the same host, the result will be equal   
tut = df_of_qseqid.loc[df_of_qseqid['level of difference'] == ""]
tut = tut.reset_index()
for d in tut['qseqid']:
    count_rows = concat_blast.loc[concat_blast['qseqid'] == d]
    count_rows = count_rows.reset_index()
    count_rows ['equal']=""
    if len(count_rows.index) == 1:
        tut ['level of difference'][tut['qseqid']== d]=  1
    else:
        for s in count_rows.index:
            if len(count_rows['spacer host taxonomy'][0]) == len(count_rows['spacer host taxonomy'][s]):
                if count_rows['spacer host taxonomy'][s] == count_rows['spacer host taxonomy'][0]:
                    count_rows ['equal'][s]=  'yes'
                else:
                   count_rows ['equal'][s]=  'no' 
            else:
                count_rows ['equal'][s]=  'no'
    if 'no' in count_rows['equal'].unique():
        tut['level of difference'][tut['qseqid']== d]=  'check'
    elif 'yes' in count_rows['equal'].unique():
        tut['level of difference'][tut['qseqid']== d]=  'equal'

print('printing with level of difference')
print(df_of_qseqid.shape)
df_of_qseqid.to_csv('plasmids_grades_results2.csv')

print('printing blanks')
print(tut.shape)
tut.to_csv('blank_results2.csv')
'''
