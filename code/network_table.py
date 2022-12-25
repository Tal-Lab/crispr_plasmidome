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
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from mobility_analysis import read_copla as mob
from mobility_analysis import mob_grades as grades

pd.set_option('display.max_columns', None)

"""
1. https://stackoverflow.com/questions/110362/how-can-i-find-the-current-os-in-python
2. https://dev.to/jakewitcher/using-env-files-for-environment-variables-in-python-applications-55a1 depending on result of 1
"""

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
blast_results = f"{resource}/BLASTp_DataBase3.zip"
mobility = f"{tables}/mobile_grades.csv"
all_ptus = f"{path}/PTUs-Mapi.xlsx"

def colors_gen(x):
    print('Generating %d colors' % x)
    all_colors = []
    i = 0
    while i < x:
        r = random.randint(0, 255)
        g = random.randint(0, 255)
        b = random.randint(0, 255)
        rgb = f'{r}, {g}, {b}'
        if rgb != '0,152,152' and rgb not in all_colors:
            all_colors.append(rgb)
            i += 1
    return all_colors

def remove_groups (taxonomy):
    try:
        if 'group' in taxonomy[5]:
            taxonomy.pop(5)
    except IndexError:
        pass
        #print(taxonomy)
    return taxonomy

def Convert (string):
    string_clean = string[1:-1].replace("'", "")
    li = list(string_clean.split(","))
    return li

def df_read():
    df = pd.read_csv(blast_results)
    df = df.loc[df['qseqid'] != df['sseqid']]
    df['spacer host taxonomy'] = df['spacer host taxonomy'].apply(Convert)
    # remove groups that are not relevent to the taxonomy hierarchy
    df['spacer host taxonomy'] = df['spacer host taxonomy'].apply(remove_groups)
    df['Superkingdom'] = df['spacer host taxonomy'].apply(lambda x: x[0] if len(x) > 0 else None)
    df['Phylum'] = df['spacer host taxonomy'].apply(lambda x: x[1] if len(x) > 1 else None)
    df['Class'] = df['spacer host taxonomy'].apply(lambda x: x[2] if len(x) > 2 else None)
    df['Order'] = df['spacer host taxonomy'].apply(lambda x: x[3] if len(x) > 3 else None)
    df['Family'] = df['spacer host taxonomy'].apply(lambda x: x[4] if len(x) > 4 else None)
    df['Genus'] = df['spacer host taxonomy'].apply(lambda x: x[5] if len(x) > 5 else None)
    df = df.rename(columns = {'spacer host species': 'Species'})
    #df = df.drop('spacer host taxonomy', axis = 1)
    return df

def df_family():
    df = df_read()
    df.Family.fillna(df['plasmid family'], inplace = True)
    df = df[['qseqid', 'Family']]
    df = df.drop_duplicates()
    return df

def df_phylum():
    df = df_read()
    df.Phylum.fillna(df['plasmid phylum'], inplace = True)
    df = df[['qseqid', 'Phylum']]
    df = df.drop_duplicates()
    return df

def color_nodes(query, level, file_name):
    pl = pd.DataFrame(query.dropna().unique(), columns=['id'])
    pl['type'] = 'plasmid'
    pl['colors'] = '0,152,152'
    pl['size'] = 40.0
    fam = pd.DataFrame(level.dropna().unique(), columns=['id'])
    fam['type'] = fam['id']
    fam['colors'] = colors_gen(len(fam))
    fam['size'] = 180.0
    df_type = pd.concat([pl,fam])
    print(df_type)
    file_name = file_name+'.csv'
    nodes_csv = f'{tables}/{file_name}'
    if not os.path.isfile(nodes_csv) or os.stat(nodes_csv).st_size == 0:
        df_type.to_csv(nodes_csv, index = False)

def network_table_family():
    df = df_family()
    df['count'] = df.groupby('Family')['Family'].transform('count')
    df.sort_values('count', ascending = False, inplace = True)
    ### setting cutoff for families with low abundance
    #df = df.loc[df['count'] >= 20]
    df = df.drop('count', axis = 1)
    df = df.reset_index(drop = True)
    print(df)
    network_csv = f'{tables}/plasmid_host_network5.csv'
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
    nodes_csv = f'{tables}/nodes5.csv'
    if not os.path.isfile(nodes_csv) or os.stat(nodes_csv).st_size == 0:
        df_type.to_csv(nodes_csv, index = False)
    return df

def network_table_phylum():
    df = df_phylum()
    df['count'] = df.groupby('qseqid')['Phylum'].transform('count')
    df.sort_values('count', ascending = False, inplace = True)
    df['edge'] = df['count'].apply(lambda x: 10 if x > 2 else 3)
    df['edge_col'] = df['edge'].apply(lambda x: '255,87,51' if x>3 else '235,235,235')
    df.loc[df['edge'] > 3, 'edge_col'] = '255,87,51'
    df.loc[df['edge'] == 3, 'edge_col'] = '235,235,235'
    print(df)
    network_csv = f'{tables}/plasmid_hostPhylum_network2.csv'
    if not os.path.isfile(network_csv) or os.stat(network_csv).st_size == 0:
        df.to_csv(network_csv, index = False)
    pl = pd.DataFrame(df['qseqid'].dropna().unique(), columns=['id'])
    pl['type'] = 'plasmid'
    pl['colors'] = '0,152,152'
    pl['size'] = 40.0
    fam = pd.DataFrame(df['Phylum'].dropna().unique(), columns=['id'])
    fam['type'] = fam['id']
    fam['colors'] = colors_gen(len(fam))
    fam['size'] = 200.0
    df_type = pd.concat([pl,fam])
    print(df_type)
    nodes_csv = f'{tables}/nodes_phylum2.csv'
    if not os.path.isfile(nodes_csv) or os.stat(nodes_csv).st_size == 0:
        df_type.to_csv(nodes_csv, index = False)
    return df

def compareFamilies(familyX, familyY):
    #print(familyX, familyY)
    #print()

    familyXname = list(familyX.keys())[0]
    familyYname = list(familyY.keys())[0]
    #print(familyXname, familyYname)

    familyXvalues = list(familyX.values())[0]
    familyYvalues = list(familyY.values())[0]

    common_list = set(familyXvalues).intersection(familyYvalues)
    #print("Common list ", common_list)
    inX = len(common_list)/len(familyXvalues)
    inY = len(common_list)/len(familyYvalues)

    result = {
        "% of common elements in family {}".format(familyXname): inX,
        "% of common elements in family {}".format(familyYname): inY,
    }

    #print(result)
    df_to_append = pd.DataFrame([[familyXname, familyYname, inX], [familyYname, familyXname, inY]], columns=['FamilyX','FamilyY','Percent'])
    #print(df_to_append)
    return df_to_append

def pivot():
    df = df_family()
    ### counting family abundance
    df['count'] = df.groupby('Family')['Family'].transform('count')
    df.sort_values('count', ascending = False, inplace = True)
    ### setting cutoff for families with low abundance
    df = df.loc[df['count'] >= 20]
    df = df.drop('count', axis = 1)
    df = df.reset_index(drop=True)
    ### grouping Families and making dictionary for pairwise comparison
    dict1 = df.groupby('Family')['qseqid'].apply(list).to_dict()
    #print(dict1)
    families = list(dict1.keys())
    #print(families)
    df_all = pd.DataFrame(columns = ['FamilyX','FamilyY','Percent'])
    for index, family in enumerate(families):
        incremetor = 1

        while index+incremetor < len(families):
            #print("\n", families[index], families[index+incremetor])

            df_index=compareFamilies(
                {families[index]: dict1[families[index]]},
                {families[index+incremetor]: dict1[families[index+incremetor]]}
            )
            df_all=df_all.append(df_index)
            incremetor += 1
    #print(df_all)
    df_pivot = pd.pivot_table(df_all, values='Percent', index='FamilyX', columns='FamilyY', fill_value = 1)
    #print(df_pivot)
    return df_pivot

def pivot_PlF():
    df = df_family()
    ### counting family abundance
    df['count'] = df.groupby('Family')['Family'].transform('count')
    df['count'] = 1
    df = df.rename(columns = {'qseqid': 'Plasmid'})
    plasmids = grades()['qseqid'].unique().tolist()  # obtaining plasmids
    all_plasmdis = df['Plasmid'].unique().tolist()
    non_reliable = [plasmid for plasmid in all_plasmdis if plasmid not in plasmids]
    print(non_reliable[0])
    plasmids.extend(non_reliable)
    df_pivot = pd.pivot_table(df, values = 'count', index = 'Plasmid', columns = 'Family', fill_value = 0)
    df_pivot = df_pivot.reindex(plasmids)
    plasmids_families = f'{tables}/plasmids_families2.xlsx'
    if not os.path.isfile(plasmids_families) or os.stat(plasmids_families).st_size == 0:
        df_pivot.to_excel(plasmids_families, index = True)

def visual():
    df = pivot()
    df = df.rename_axis("Families")
    df = df.rename_axis("Families", axis = "columns")
    ### creating clustermap
    cg=sns.clustermap(df,
                      metric = "euclidian",
                      method = "ward",
                      xticklabels = True,
                      yticklabels = True,
                      cmap = 'Blues',
                      cbar_pos=(0, .3, .03, .4),
                      cbar_kws={'ticks':[0,1],
                                'label': 'Plasmid Percentage'}
                      )
    cg.ax_row_dendrogram.set_visible(False)
    cg.ax_col_dendrogram.set_visible(False)
    svg_name = 'Clustermap.svg'
    svg_file = f'{visuals}/{svg_name}'
    png_name = 'Clustermap.png'
    png_file = f'{visuals}/{png_name}'
    if not os.path.isfile(svg_file) and not os.path.isfile(png_file):
        plt.savefig(svg_file, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
        plt.savefig(png_file, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    plt.show()
    ### creating heatmap
    fig = sns.heatmap(df,
                      xticklabels = True,
                      yticklabels = True,
                      cmap = 'Blues',
                      cbar_kws={'ticks':[0,1],
                                'label': 'Plasmid Percentage'}
                      )
    svg_name = 'Heatmap.svg'
    svg_file = f'{visuals}/{svg_name}'
    png_name = 'Heatmap.png'
    png_file = f'{visuals}/{png_name}'
    if not os.path.isfile(svg_file) and not os.path.isfile(png_file):
        plt.savefig(svg_file, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
        plt.savefig(png_file, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    plt.show()

def mob_grades():
    mob_df = mob()
    mob_df = mob_df.groupby('Plasmid').agg(MOB = ('MOB', ','.join)).reset_index()
    grade_df = grades()
    return [mob_df, grade_df]

def df_family_top30():
    df = df_read()
    df = df[['qseqid', 'Family', 'Phylum']]
    df = df.drop_duplicates()
    ### getting top 30 abundant families and merging others into OTHER
    df['count'] = df.groupby('Family')['Family'].transform('count')
    df.sort_values('count', ascending = False, inplace = True)
    top30_fams = df['Family'].unique().tolist()[:30]
    df['Family'] = df['Family'].apply(lambda x: x if x in top30_fams else 'Other')
    ### counting number of conntections for each plasmid
    df = df.drop('count', axis = 1)
    df['count'] = df.groupby('qseqid')['Family'].transform('count')
    df.sort_values('count', ascending = False, inplace = True)
    df['edge'] = df['count'].apply(lambda x: x*2)
    df = df.reset_index(drop = True)
    print(df)
    ### adding mobility info for plasmids
    df_withMOB = df.merge(mob_grades()[0], how = 'left', left_on = 'qseqid', right_on = 'Plasmid')
    df_withMOB.drop(['Plasmid'], axis=1, inplace = True)
    df_withMOB['MOB'] = df_withMOB['MOB'].fillna('MOB-')

    ### adding host-range grade info for plasmids
    grades_df = mob_grades()[1]
    grades_df = grades_df[['qseqid', 'level of difference']]
    df_withMOB = df_withMOB.merge(grades_df, how = 'left', on = 'qseqid')
    #print(df_withMOB)
    network_csv = f'{tables}/plasmid_host_network_top30_3.csv'
    if not os.path.isfile(network_csv) or os.stat(network_csv).st_size == 0:
        df_withMOB.to_csv(network_csv, index = False)
    pl = pd.DataFrame(df['qseqid'].dropna().unique(), columns=['id'])
    pl['type'] = 'plasmid'
    pl['colors'] = '0,152,152'
    pl['size'] = 40.0
    fam = pd.DataFrame(df['Family'].dropna().unique(), columns=['id'])
    fam['type'] = fam['id']
    fam['colors'] = colors_gen(len(fam))
    fam['size'] = 180.0
    df_type = pd.concat([pl,fam])
    print(df_type)
    nodes_csv = f'{tables}/nodes_top30_3.csv'
    if not os.path.isfile(nodes_csv) or os.stat(nodes_csv).st_size == 0:
        df_type.to_csv(nodes_csv, index = False)
    return df

def df_phylum2():
    df = df_phylum()
    df = df.drop_duplicates()
    df['count'] = df.groupby('qseqid')['Phylum'].transform('count')
    df.sort_values('count', ascending = False, inplace = True)
    df['edge'] = df['count'].apply(lambda x: 10 if x > 2 else 3)
    df['edge_col'] = df['edge'].apply(lambda x: '255,87,51' if x>3 else '235,235,235')
    df.loc[df['edge'] > 3, 'edge_col'] = '255,87,51'
    df.loc[df['edge'] == 3, 'edge_col'] = '235,235,235'
    print(df)
    ### adding mobility info for plasmids
    df_withMOB = df.merge(mob_grades()[0], how = 'left', left_on = 'qseqid', right_on = 'Plasmid')
    df_withMOB.drop(['Plasmid'], axis = 1, inplace = True)
    df_withMOB['MOB'] = df_withMOB['MOB'].fillna('MOB-')

    ### adding host-range grade info for plasmids
    grades_df = mob_grades()[1]
    grades_df = grades_df.rename(columns={'MOB':'Mobility'})
    df_withMOB = df_withMOB.merge(grades_df, how = 'left', on = 'qseqid')
    print(df_withMOB)
    network_csv = f'{tables}/plasmid_hostPhylum_network3.csv'
    if not os.path.isfile(network_csv) or os.stat(network_csv).st_size == 0:
        df_withMOB.to_csv(network_csv, index = False)
    pl = pd.DataFrame(df['qseqid'].dropna().unique(), columns=['id'])
    pl['type'] = 'plasmid'
    pl['colors'] = '0,152,152'
    pl['size'] = 40.0
    fam = pd.DataFrame(df['Phylum'].dropna().unique(), columns=['id'])
    fam['type'] = fam['id']
    fam['colors'] = colors_gen(len(fam))
    fam['size'] = 200.0
    df_type = pd.concat([pl,fam])
    print(df_type)
    nodes_csv = f'{tables}/nodes_phylum3.csv'
    if not os.path.isfile(nodes_csv) or os.stat(nodes_csv).st_size == 0:
        df_type.to_csv(nodes_csv, index = False)
    return df

def ptus ():
    df = pd.read_excel(all_ptus, header = 1)
    df = df[["AccessionVersion", "PTU", "Hrange (1)"]]
    df = df.rename(columns = {'AccessionVersion': 'qseqid',"Hrange (1)":"Hrange" })
    #plasmids_ptu = df['AccessionVersion'].unique().tolist()
    df_fam = df_family()
    df_fam = df_fam.merge(df, how = 'left', on = 'qseqid')
    df_no_nan = df_fam.loc[~df_fam['PTU'].isnull()]
    df_no_nan = df_no_nan.loc[df_fam['PTU'] != '-']
    df_fam.PTU.fillna('missing', inplace = True)
    #print(df_no_nan)
    #print(df_fam['PTU'].unique())
    return df_fam, df_no_nan

def ptu_network(dataset):
    if dataset == 'no-nan':
        df = ptus()[1]
    elif dataset == 'miss':
        df = ptus()[0]
    else:
        print("specified dataset is not valid")

    ### counting number of conntections for each PTU
    df = df[['Family','PTU']].drop_duplicates()
    df['count'] = df.groupby('PTU')['Family'].transform('count')
    df.sort_values('count', ascending = False, inplace = True)
    #df['edge'] = 10
    df = df.reset_index(drop = True)
    network_csv = f'{tables}/PTU_host_network.csv'
    if not os.path.isfile(network_csv) or os.stat(network_csv).st_size == 0:
        df.to_csv(network_csv, index = False)
    color_nodes(df['PTU'], df['Family'], 'nodes_PTU_Fam2')
    return df

def PTU_family_phylum():
    df = df_read()
    df.Family.fillna(df['plasmid family'], inplace = True)
    df.Phylum.fillna(df['plasmid phylum'], inplace = True)
    df.Class.fillna(df['plasmid class'], inplace = True)
    df.Order.fillna(df['plasmid order'], inplace = True)
    df = df[['Family', 'Order', 'Class', 'Phylum']]
    df = df.drop_duplicates()
    df_ptu= ptu_network('no-nan')
    df_ptu_family = df_ptu.merge(df, how = 'left', on ='Family')
    network_csv = f'{tables}/PTU_hostLin_network.csv'
    if not os.path.isfile(network_csv) or os.stat(network_csv).st_size == 0:
        df_ptu_family.to_csv(network_csv, index = False)

def ptu_grades():
    df = ptus()[1]
    df_blast = df_read()

    df_ptus = df.merge(df_blast, how = 'left', on = 'qseqid')
    list_of_ptus = df_ptus['PTU'].unique()

    # creating new dataframe
    df_of_ptu = pd.DataFrame(df_ptus['PTU'].unique())
    df_of_ptu['level of difference'] = ""
    df_of_ptu.rename(columns = {0: 'PTU'}, inplace = True)
    # go over qseqid list and create a dataframe from all result of the specific id in blast
    for i in list_of_ptus:
        print(i)
        rslt_df = df_ptus.loc[df_ptus['PTU'] == i]
        rslt_df = rslt_df.reset_index()
        rslt_df['level of difference'] = ""
        # go over rows of the temporary datafrme of the specific id and check if there are results for different hosts(taxonomy), if so check the level of difference
        for e in rslt_df.index:
            try:
                print('Trying for %s in %s' % (e, i))
                if len(rslt_df['spacer host taxonomy']) != len(rslt_df['spacer host taxonomy'][e]):
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
                    elif rslt_df['Species'][0] != rslt_df['Species'][e]:
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
                        elif rslt_df['Species'][0] != rslt_df['Species'][e]:
                            rslt_df.loc[e, 'level of difference'] = 'Species'
                    else:
                        if rslt_df['Species'][0] != rslt_df['Species'][e]:
                            rslt_df.loc[e, 'level of difference'] = 'Species'
            except IndexError:
                if rslt_df['Species'][0] != rslt_df['Species'][e]:
                    rslt_df.loc[e, 'level of difference'] = 'Species'
                    # check the most far taxonomy difference
        if 'Phylum' in rslt_df['level of difference'].unique():
            df_of_ptu['level of difference'][df_of_ptu['PTU'] == i] = "6"
        elif 'Class' in rslt_df['level of difference'].unique():
            df_of_ptu['level of difference'][df_of_ptu['PTU'] == i] = "5"
        elif "Order" in rslt_df['level of difference'].unique():
            df_of_ptu['level of difference'][df_of_ptu['PTU'] == i] = "4"
        elif "Family" in rslt_df['level of difference'].unique():
            df_of_ptu['level of difference'][df_of_ptu['PTU'] == i] = "3"
        elif "Genus" in rslt_df['level of difference'].unique():
            df_of_ptu['level of difference'][df_of_ptu['PTU'] == i] = "2"
        elif "Species" in rslt_df['level of difference'].unique():
            df_of_ptu['level of difference'][df_of_ptu['PTU'] == i] = "1"
    tut = df_of_ptu[df_of_ptu['level of difference'] == ""]
    tut = tut.reset_index()
    for d in tut['PTU']:
        count_rows = df_of_ptu.loc[df_of_ptu['PTU'] == d]
        count_rows = count_rows.reset_index()
        count_rows['equal'] = ""
        if len(count_rows.index) == 1:
            tut['level of difference'][tut['PTU'] == d] = 1
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
            tut['level of difference'][tut['PTU'] == d] = 'check'
        elif 'yes' in count_rows['equal'].unique():
            tut['level of difference'][tut['PTU'] == d] = 'equal'
        df_1 = tut.loc[tut['level of difference'] != 'equal']
        plasmids_to_remove = df_1['PTU'].unique().tolist()
        # print(plasmids_to_remove)
        df_of_ptu['PTU'] = df_of_ptu['PTU'].apply(lambda x: x if x not in plasmids_to_remove else 'remove')
        df_of_ptu = df_of_ptu.loc[df_of_ptu['PTU'] != 'remove']
        df_of_ptu['level of difference'].fillna(1, inplace = True)

    PTU_hostrange = f'{tables}/PTU_HostRange.csv'
    if not os.path.isfile(PTU_hostrange) or os.stat(PTU_hostrange).st_size == 0:
        df_of_ptu.to_csv(PTU_hostrange, index = True)


#ptu_grades()
def pivot_PTUs():
    #df = df_family()
    df = ptus()[0]
    df = df[['qseqid', 'Family','PTU']]
    df.loc[df['PTU']=='missing', 'PTU'] = '-'
    print(df)
    ### counting family abundance
    dict1 = df.groupby('Family')['PTU'].apply(list).to_dict()
    # print(dict1)
    families = list(dict1.keys())
    # print(families)
    df_all = pd.DataFrame(columns = ['FamilyX', 'FamilyY', 'Percent'])
    for index, family in enumerate(families):
        incremetor = 1

        while index + incremetor < len(families):
            # print("\n", families[index], families[index+incremetor])

            df_index = compareFamilies(
                {families[index]: dict1[families[index]]},
                {families[index + incremetor]: dict1[families[index + incremetor]]}
            )
            df_all = df_all.append(df_index)
            incremetor += 1
    #df['count'] = 1
    df_all = df_all.rename(columns = {'FamilyX': 'Family'})
    df_merged = df_all.merge(df[['Family', 'PTU']],how = 'left', on = 'Family')
    print(df_merged)
    df_pivot = pd.pivot_table(df_merged, values='Percent', index=['PTU','Family'], columns='FamilyY', fill_value = 1)
    print(df_pivot)
    PTU_pivot_families = f'{tables}/PTU_families2.xlsx'
    if not os.path.isfile(PTU_pivot_families) or os.stat(PTU_pivot_families).st_size == 0:
        df_pivot.to_excel(PTU_pivot_families, index = True)
    PTU_families = f'{tables}/PTU_families_notpivot2.csv'
    if not os.path.isfile(PTU_families) or os.stat(PTU_families).st_size == 0:
        df_merged.to_csv(PTU_families, index = True)

#ptu_network('no-nan')
#pivot_PTUs()
#PTU_family_phylum()
#df_family_top30()
#df_phylum2()
#pivot_PlF()
#visual()

