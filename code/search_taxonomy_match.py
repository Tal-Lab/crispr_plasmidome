# -*- coding: utf-8 -*-
"""
Created on Sun Dec 11 19:02:23 2022

@author: sivan
"""

import pandas as pd

plasmid_blast = pd.read_csv (r"C:\sivan\לימודים\שנה ג\פרויקט מחקר בפייתון\BLASTp_DataBase3.csv", sep= ",",  header = 0)
#remove rows where qseqid == sseqid
plasmid_blast = plasmid_blast[plasmid_blast.qseqid != plasmid_blast.sseqid]

#convert string to list
def Convert(string):
    string_clean = string [1:-1].replace("'", "")    
    li = list(string_clean.split(","))
    return li

#apply conver function on 'spacer host taxonomy' column
plasmid_blast['spacer host taxonomy'] = plasmid_blast['spacer host taxonomy'].apply(Convert)




#remove groups that are not relevent to the taxonomy hierarchy 
def remove_groups(taxonomy):
    try:
        if 'group' in taxonomy[5]:
            taxonomy.pop(5)
    except IndexError:
        print(taxonomy)
    return taxonomy


#apply remove_groups function on 'spacer host taxonomy' column
plasmid_blast['spacer host taxonomy'] = plasmid_blast['spacer host taxonomy'].apply(remove_groups)


plasmid_blast.to_csv("clean_up.csv")
plasmid_blast = pd.read_csv (r"C:\sivan\לימודים\clean_up.csv", sep= ",",  header = 0)

#clean up column before comparison
def remove_unnecessary_charecters (cell):
    return cell[2:-1]

#apply remove_unnecessary_charecters function on spacer superkingdom column
plasmid_blast['spacer superkingdom'] = plasmid_blast['spacer superkingdom'].apply(remove_unnecessary_charecters)

#clean up column before comparison
def remove_space_and_other (string):
    try:
        string = string.replace(' ' , '')
    except AttributeError:
        return np.nan
        
    return string.replace(']' , '')

#apply remove_space_and_other function on columns before comparison
plasmid_blast['spacer genus'] = plasmid_blast['spacer genus'].apply(remove_space_and_other)
plasmid_blast['spacer family'] = plasmid_blast['spacer family'].apply(remove_space_and_other)
plasmid_blast['spacer order'] = plasmid_blast['spacer order'].apply(remove_space_and_other)
plasmid_blast['spacer class'] = plasmid_blast['spacer class'].apply(remove_space_and_other)
plasmid_blast['spacer phylum'] = plasmid_blast['spacer phylum'].apply(remove_space_and_other)


#search in every hit, the nearest alignment to known host
for i in plasmid_blast.index:
    if plasmid_blast['spacer host species'][i] == plasmid_blast['plasmid species'][i]:
        plasmid_blast['match'][i] = 'Species'
    elif plasmid_blast['spacer genus'][i] == plasmid_blast['plasmid genus'][i]:
        plasmid_blast['match'][i] = 'Genus'
    elif plasmid_blast['spacer family'][i] == plasmid_blast['plasmid family'][i]:
        plasmid_blast['match'][i] = 'Family'
    elif plasmid_blast['spacer order'][i] == plasmid_blast['plasmid order'][i]:
        plasmid_blast['match'][i] = 'Order'
    elif plasmid_blast['spacer class'][i] == plasmid_blast['plasmid class'][i]:
        plasmid_blast['match'][i] = 'Class'
    elif plasmid_blast['spacer phylum'][i] == plasmid_blast['plasmid phylum'][i]:
        plasmid_blast['match'][i] = 'Phylum'
    elif plasmid_blast['spacer superkingdom'][i] == plasmid_blast['plasmid superkingdom'][i]:
        plasmid_blast['match'][i] = 'Superkingdom'

plasmid_blast.to_csv('match_update.csv')
