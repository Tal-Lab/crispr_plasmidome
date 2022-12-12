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

def df_family():
    df = pd.read_csv(blast_results)
    df = df[['qseqid', 'sseqid', 'ratio',  'plasmid family','spacer host taxonomy']]
    df = df.loc[df['qseqid'] != df['sseqid']]
    df['spacer host taxonomy'] = df['spacer host taxonomy'].apply(lambda x: ast.literal_eval(x))
    df[['Kingdom','Phylum','Class','Order','Family', 'Genus','else1', 'else2']] = pd.DataFrame(df['spacer host taxonomy'].tolist())
    df.Family.fillna(df['plasmid family'], inplace = True)
    df = df[['qseqid', 'Family']]
    df = df.drop_duplicates()
    df['count'] = df.groupby('Family')['Family'].transform('count')
    df.sort_values('count', ascending = False, inplace = True)
    ### setting cutoff for families with low abundance
    df = df.loc[df['count'] >= 20]
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
#df_family()
def df_phylum():
    df = pd.read_csv(blast_results)
    df = df[['qseqid', 'sseqid', 'ratio',  'plasmid phylum','spacer host taxonomy']]
    df = df.loc[df['qseqid'] != df['sseqid']]
    df['spacer host taxonomy'] = df['spacer host taxonomy'].apply(lambda x: ast.literal_eval(x))
    df[['Kingdom','Phylum','Class','Order','Family', 'Genus','else1', 'else2']] = pd.DataFrame(df['spacer host taxonomy'].tolist())
    df.Phylum.fillna(df['plasmid phylum'], inplace = True)
    df = df[['qseqid', 'Phylum']]
    df = df.drop_duplicates()
    print(df)
    network_csv = f'{tables}/plasmid_hostPhylum_network.csv'
    if not os.path.isfile(network_csv) or os.stat(network_csv).st_size == 0:
        df.to_csv(network_csv, index = False)
    pl = pd.DataFrame(df['qseqid'].dropna().unique(), columns=['id'])
    pl['type'] = 'plasmid'
    pl['colors'] = '0,152,152'
    pl['size'] = 20.0
    fam = pd.DataFrame(df['Phylum'].dropna().unique(), columns=['id'])
    fam['type'] = fam['id']
    fam['colors'] = colors_gen(len(fam))
    fam['size'] = 150.0
    df_type = pd.concat([pl,fam])
    print(df_type)
    nodes_csv = f'{tables}/nodes_phylum.csv'
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

#visual()