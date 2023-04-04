"""
Created on 04/04/2023

Author: Lucy Androsiuk
"""

### importing modules
import os, re
import pandas as pd
from Bio import Entrez
from Bio import SeqIO
from pathlib import Path
import logging

#pd.set_option('display.max_colwidth', None)
### paths
email = "androsiu@post.bgu.ac.il"  # Tell NCBI who you are
# uncomment relevant path to OS
# Windows
#path = r"C:\Users\Lucy\iCloudDrive\Documents\bengurion\Project students\Sivan_project"
# macOS
#path = r"/Users/lucyandrosiuk/Documents/bengurion/Project students/Sivan_project"
# Cluster
path = r"/gpfs0/tals/projects/Analysis/Lucy_plasmidome/Plasmidome/CRISPR"

# working directories
visuals = f"{path}/visualisations"
tables = f"{path}/data_calculations"
resource = f"{path}/res"
Path(visuals).mkdir(parents=True, exist_ok=True)

# working files
plasmids_df = f"{path}/match_update_90.csv"
all_features_csv = f'{tables}/all_features.csv'

logging.basicConfig(filename=f'{path}/log_all_features.log', format = '%(asctime)s - %(levelname)s - %(message)s', level=logging.DEBUG)
logging.debug("I'm inside the python file")

def get_feature_start_end(feature):
    location = feature.location
    strand = feature.strand

    sublocations = location.parts
    num_sublocations = len(sublocations)

    if num_sublocations == 1:
        start = sublocations[0].start
        end = sublocations[0].end
    else:
        if strand == 1:
            start = sublocations[0].start
            end = sublocations[-1].end
        else:
            start = sublocations[-1].start
            end = sublocations[0].end
    return start, end

def search_Entrez(x):
    logging.debug('Entered search_Entrez for plasmid %s' % x)
    Entrez.email = email
    Entrez.sleep_between_tries = 20  # Tell NCBI delay, in seconds, before retrying a request
    try:
        # retrieve the records from NCBI in GenBank format
        handle = Entrez.efetch(db = "nuccore", id = x, rettype = "gb", retmode = "text")
        record = SeqIO.read(handle, "gb")
        # setting lists for future df
        pl_ids = []
        ids = []
        p_types = []
        p_names = []
        seqs = []
        p_starts = []
        p_ends = []
        go_funcs = []

        # extract the features
        for feature in record.features:
            if feature.type == 'gene' or feature.type == 'source':
                pass
            elif feature.type == 'CDS' or feature.type == 'mobile_element':
                pl_ids.append(x)
                #### feature type
                p_type = feature.type
                p_types.append(p_type)
                #### location
                start, end = get_feature_start_end(feature)
                p_starts.append(start)
                p_ends.append(end)

                #### qualifiers

                ## feature_id
                feature_id = feature.qualifiers.get("protein_id", [""])[0]
                mobile_element_value = feature.qualifiers.get("mobile_element_type", [""])[0]
                note = feature.qualifiers.get("note", [""])[0]

                if feature_id:
                    # print('ID:', feature_id)
                    feature_id = feature_id
                elif mobile_element_value:
                    # print("ID:", 'mobile_element_%d' % i)
                    feature_id = note
                else:
                    # print("ID:", 'pseudogene_%d' % i)
                    feature_id = note
                ids.append(feature_id)

                ## product name
                product_value = feature.qualifiers.get("product", [""])[0]
                pseudogene_value = feature.qualifiers.get("pseudogene", [""])[0]
                if product_value:
                    # print("Product:", product_value)
                    p_name = product_value
                elif mobile_element_value:
                    # print("Product:", mobile_element_value)
                    p_name = mobile_element_value
                elif pseudogene_value:
                    # print("Product:", pseudogene_value)
                    p_name = pseudogene_value
                else:
                    p_name = note
                p_names.append(p_name)

                ## GO function
                Go_function_value = feature.qualifiers.get("Go_function", [""])[0]
                if Go_function_value:
                    go_func = Go_function_value
                    # print("Go Function:", Go_function_value)
                else:
                    # print("Go Function:", 'missing')
                    go_func = 'missing'
                go_funcs.append(go_func)

                ## sequence
                sequence = feature.qualifiers.get("translation", [""])[0]
                if sequence:
                    sequence = sequence
                else:
                    sequence = '-'

                seqs.append(sequence)

            else:
                pl_ids.append(x)
                #### feature type
                p_type = feature.type
                p_types.append(p_type)
                #### location
                start, end = get_feature_start_end(feature )
                p_starts.append(start)
                p_ends.append(end)
                #### qualifiers

                ## feature_id
                feature_id = feature.qualifiers.get("protein_id", [""])[0]
                mobile_element_value = feature.qualifiers.get("mobile_element_type", [""])[0]
                note = feature.qualifiers.get("note", [""])[0]

                if feature_id:
                    feature_id = feature_id
                else:
                    feature_id = note
                ids.append(feature_id)

                ## product name
                product_value = feature.qualifiers.get("product", [""])[0]
                pseudogene_value = feature.qualifiers.get("pseudogene", [""])[0]

                if product_value:
                    # print("Product:", product_value)
                    p_name = product_value
                elif mobile_element_value:
                    # print("Product:", mobile_element_value)
                    p_name = mobile_element_value
                elif pseudogene_value:
                    # print("Product:", pseudogene_value)
                    p_name = pseudogene_value
                else:
                    p_name = note
                p_names.append(p_name)

                ## GO function
                Go_function_value = feature.qualifiers.get("Go_function", [""])[0]
                if Go_function_value:
                    go_func = Go_function_value
                    # print("Go Function:", Go_function_value)
                else:
                    # print("Go Function:", 'missing')
                    go_func = 'missing'
                go_funcs.append(go_func)

                ## sequence
                sequence = feature.qualifiers.get("translation", [""])[0]
                if sequence:
                    sequence = sequence

                else:
                    sequence = '-'

                seqs.append(sequence)
        df = pd.DataFrame({'qseqid': pl_ids, 'Type': p_types, 'ID': ids, 'P_Name': p_names, 'GO_function': go_funcs,
                           'P_start': p_starts, 'P_end': p_ends, 'Sequence': seqs})
        logging.debug("############## Features for %s obtained ##############", x)
        if not os.path.isfile(all_features_csv) or os.stat(all_features_csv).st_size == 0:
            df.to_csv(all_features_csv, mode = 'a', index = True, header = True)
        else:
            df.to_csv(all_features_csv, mode = 'a', index = True, header = False)
        # print(record.annotations["structured_comment"])
        handle.close()
    except Exception as e:
        logging.debug(e)

df_range = pd.read_csv(plasmids_df, sep = ',', header = 0, index_col = 0)
plasmids_list = df_range['qseqid'].unique().tolist()
for plasmid in plasmids_list:
    search_Entrez(plasmid)