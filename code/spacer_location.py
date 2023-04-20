"""
Created on 30/03/2023

Author: Lucy Androsiuk
"""

### importing modules
import pandas as pd
import numpy as np
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os, re
from pathlib import Path
import logging
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
import seaborn as sns
from scipy import stats

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

#pd.set_option('display.max_colwidth', None)
### paths
email = "androsiu@post.bgu.ac.il"  # Tell NCBI who you are
# uncomment relevant path to OS
# Windows
#path = r"C:\Users\Lucy\iCloudDrive\Documents\bengurion\Project students\Sivan_project"
# macOS
path = r"/Users/lucyandrosiuk/Documents/bengurion/Project students/Sivan_project"
# Cluster
#path = r"/gpfs0/tals/projects/Analysis/Lucy_plasmidome/Plasmidome/CRISPR"

# working directories
visuals = f"{path}/visualisations"
tables = f"{path}/data_calculations"
resource = f"{path}/res"
Path(visuals).mkdir(parents=True, exist_ok=True)
features = f'{tables}/all_features.csv'

def work_files(cutoff):
    #host_grades = f"{path}/cutoffs/id_{cutoff}/plasmids_grades_results_{cutoff}.csv"
    #blank = f"{path}/cutoffs/id_{cutoff}/blank_results_{cutoff}.csv"
    match_update = f"{path}/match_update_{cutoff}.csv"
    all_info = f"{tables}/all_info_{cutoff}.csv"
    #return host_grades, blank, match_update, all_info
    return match_update, all_info

def spacers_locations(cutoff):
    df_spacer = pd.read_csv(work_files(cutoff)[0], sep = ',', header = 0, index_col = 0)
    df_spacer = df_spacer[['qseqid', 'sseqid', 'length', 'ratio', 'qstart', 'qend', 'match']]

    df_features = pd.read_csv(features, sep = ',', header = 0, index_col = 0)
    list_of_plasmids = df_features['qseqid'].unique().tolist()

    all_plasmids = df_range['qseqid'].unique().tolist()
    difference_result = []
    for pl in all_plasmids:
        if pl not in list_of_plasmids:
            difference_result.append(pl)

    '''
    with open(f'{tables}/plasmids_missing.txt', 'w') as fp:
        fp.write('\n'.join(difference_result))
    '''
    # create an empty list to store the matches
    all_matches = []

    # loop over each genome_id
    for genome_id in list_of_plasmids:
        print(genome_id)
        # get the features and matches for the current genome_id
        features_df = df_features[df_features['qseqid'] == genome_id]
        # print(features_df)
        matches_df = df_spacer[df_spacer['qseqid'] == genome_id]
        # print(matches_df)

        # loop over each feature row
        for feature_index, feature_row in features_df.iterrows():
            feature_start = feature_row['P_start']
            feature_end = feature_row['P_end']

            # loop over each match row
            for match_index, match_row in matches_df.iterrows():
                match_start = match_row['qstart']
                match_end = match_row['qend']
                # check if the match falls within the feature coordinates
                if (match_start >= feature_start) and (match_end <= feature_end):
                    # add the match to the list of matches that fall within features
                    all_matches.append({
                        'genome_id': genome_id,
                        'Feature_Type': feature_row['Type'],
                        'Feature_ID': feature_row['ID'],
                        'Feature_Name': feature_row['P_Name'],
                        'GO_function': feature_row['GO_function'],
                        'Spacer_ID': match_row['sseqid'],
                        'Match_Start': match_row['qstart'],
                        'Match_End': match_row['qend'],
                        'Within_Feature': True})
                else:
                    # add the match to the list of all matches, even if it doesn't fall within proteins or pseudo_proteins
                    all_matches.append({
                        'genome_id': genome_id,
                        'Feature_Type': '-',
                        'Feature_ID': '-',
                        'Feature_Name': '-',
                        'GO_function': '-',
                        'Spacer_ID': match_row['sseqid'],
                        'Match_Start': match_row['qstart'],
                        'Match_End': match_row['qend'],
                        'Within_Feature': False})
    print(all_matches)

    # create a new dataframe with all matches
    df_all_matches = pd.DataFrame(all_matches)
    print(df_all_matches.sample(n=200))

    # write the results to a new CSV file
    df_all_matches.to_csv(f'{tables}/all_matchesF_{cutoff}.csv', index = True)


def spacers_loc_HR(cutoff):
    df= pd.read_csv(f'{tables}/all_matchesF_{cutoff}.csv', header = 0, index_col = 0)
    df_range = pd.read_csv(work_files(cutoff)[1], sep = ',', header = 0, index_col = 0)
    df_range = df_range[['qseqid', 'level of difference', 'MOB']]
    print(df.columns)
    df_merged = df.merge(df_range, left_on = 'genome_id', right_on='qseqid')
    df_merged = df_merged.drop('qseqid', axis=1)
    print(df_merged.sample(n=50))
    df_merged.to_csv(f'{tables}/all_matches_HR_{cutoff}.csv', index = True)

def read_loc(cutoff):
    df = pd.read_csv(f'{tables}/all_matches_HR_{cutoff}.csv', header = 0, index_col = 0)
    df = df.drop_duplicates()
    df = df.loc[df['genome_id'] != df['Spacer_ID']]
    print(df.columns)
    #print(df.sample(100))
    df_true = df.loc[df['Within_Feature'] == True]
    df_repeat = df.loc[df['Feature_Type']=='repeat_region']
    df_mobile = df_true.loc[df_true['Feature_Type']=='mobile_element']
    #df_mobile.to_csv(f'{tables}/mobile_element_{cutoff}.csv', index = True)
    #print(df_repeat.sample(100))
    #print(df_true['Feature_Name'].unique())
    print(df_true['Feature_Name'].value_counts(normalize = True))
    print(df_true['Feature_Type'].value_counts())
    feature_frequency = df_true['Feature_Name'].value_counts(normalize = True)
    #feature_frequency.to_csv(f'{tables}/feature_freq_{cutoff}.csv', index = True)
    print(len(df_true))
    print(len(df))
    return df

def get_dupl(cutoff):
    df = read_loc(cutoff)
    df = df.drop_duplicates()
    dupl = df[df[['genome_id', 'Spacer_ID']].duplicated()]
    dupl.to_csv(f'{tables}/Multiple_{cutoff}.csv', index = True)


def get_stats(cutoff):
    df = pd.read_csv(work_files(cutoff)[0], header = 0, index_col = 0)
    print(df.head)
    df = df[['qseqid','sseqid', 'qstart', 'qend']]
    df = df.drop_duplicates()
    dupl = df[df[['qseqid','sseqid']].duplicated()]
    print(dupl['qseqid'].unique())

    ### counting statistcs for hits per plasmid
    df_hits_plasmid = df['qseqid'].value_counts().reset_index(name='Number of hits')
    df_hits_plasmid.to_csv(f'{tables}/Hits_per_Plasmid_{cutoff}.csv', index = True)
    # print the statistics
    print(df_hits_plasmid.describe())

    # plot the distribution of hits
    sns.histplot(data = df_hits_plasmid, bins = 30, x = 'Number of hits')

    # set the title and x-axis label
    plt.title('Distribution of hits for each plasmid')
    plt.xlabel('Number of hits')
    svg_name = f"HitsPerPlasmidHisto_{cutoff}.svg"
    svg_dir = f'{visuals}/{svg_name}'
    png_name = f"HitsPerPlasmidHisto_{cutoff}.png"
    png_dir = f'{visuals}/{png_name}'
    # plt.autoscale()
    plt.savefig(svg_dir, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
    plt.savefig(png_dir, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    # display the plot
    plt.show()
    ### counting statistcs for spacers per plasmid
    df_spacers_plasmid = df.groupby('qseqid')['sseqid'].nunique().reset_index(name='Number of spacers')
    df_spacers_plasmid.to_csv(f'{tables}/Spacers_per_Plasmid_{cutoff}.csv', index = True)
    print(df_spacers_plasmid.describe())
    # plot the distribution of hits
    # plot the statistics
    sns.histplot(data = df_spacers_plasmid, bins = 30, x='Number of spacers')
    # set the title and x-axis label
    plt.title('Distribution of number of spacers for each plasmid')
    plt.xlabel('Number of spacers')
    svg_name = f"SpacersPerPlasmidHisto_{cutoff}.svg"
    svg_dir = f'{visuals}/{svg_name}'
    png_name = f"SpacersPerPlasmidHisto_{cutoff}.png"
    png_dir = f'{visuals}/{png_name}'
    # plt.autoscale()
    plt.savefig(svg_dir, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
    plt.savefig(png_dir, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    # display the plot
    plt.show()
    ### counting statistcs for plasmids per spacer
    df_plasmid_spacer = df.groupby('sseqid')['qseqid'].nunique().reset_index(name='Number of plasmids')
    df_plasmid_spacer.to_csv(f'{tables}/Plasmids_Per_Spacer_{cutoff}.csv', index = True)
    print(df_plasmid_spacer.describe())
    # plot the distribution of hits
    sns.histplot(data = df_plasmid_spacer, bins =  30, x= 'Number of plasmids')
    # set the title and x-axis label
    plt.title('Distribution of plasmids for each spacer')
    plt.xlabel('Number of plasmids')
    svg_name = f"PlasmidsPerSpacerHisto_{cutoff}.svg"
    svg_dir = f'{visuals}/{svg_name}'
    png_name = f"PlasmidsPerSpacerHisto_{cutoff}.png"
    png_dir = f'{visuals}/{png_name}'
    # plt.autoscale()
    plt.savefig(svg_dir, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
    plt.savefig(png_dir, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    # display the plot
    plt.show()

def cutoff():
    cutoff = [90, 95, 100]
    for i in cutoff:
        get_stats(i)

#cutoff()
