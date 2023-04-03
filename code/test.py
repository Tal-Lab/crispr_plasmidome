# -*- coding: utf-8 -*-
"""
Created on 02/04/2023 15:04

Author: Lucy
"""

import pandas as pd
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import ExactPosition, FeatureLocation, CompoundLocation
import os, re


email = "androsiu@post.bgu.ac.il"  # Tell NCBI who you are

plasmids = ['AP019704.1', 'NZ_LT719075.1', 'NZ_CP050379.1']

def search_Entrez (x):
    '''getting orfeome for each spacer-barer id'''
    print("Entering search_Entrez function with plasmid ", x)
    Entrez.email = email  # Tell NCBI who you are
    Entrez.sleep_between_tries = 20  #Tell NCBI delay, in seconds, before retrying a request
    try:
        print("############## Obtaining proteins for %s ##############", x)
        handle = Entrez.efetch(db="nuccore", id= x, rettype="fasta_cds_aa", retmode="text")
        records = SeqIO.parse(handle, "fasta")
        ids = []
        p_names = []
        seqs = []
        p_starts = []
        p_ends = []
        for record in records:
            if record.description.__contains__('[pseudo=true]'):
                print("!!!!!!!!!!!!!! PSEUDO")
                #print(record)
                id = record.id[4:]
                ids.append(id)
                description = record.description
                parts = description.split('[')
                #print(parts)
                p_loc_info = parts[-2]
                p_loc = re.findall(r'\d+', p_loc_info)
                p_start = int(p_loc[0])
                p_end = int(p_loc[-1])
                p_starts.append(p_start)
                p_ends.append(p_end)
                p_name = [match[8:-2] for match in parts if 'protein' in match]
                if p_name == []:
                    p_name = [match[11:-2] for match in parts if 'pseudogene' in match]
                print(p_name)
                p_names.append(p_name)
                sequence = str(record.seq)
                seqs.append(sequence)
            else:
                pass

            ''' 
            
            description = record.description
            parts = description.split('[')
            p_loc_info = parts[-2]
            p_loc = re.findall(r'\d+', p_loc_info)
            p_start = int(p_loc[0])
            p_end = int(p_loc[-1])
            p_starts.append(p_start)
            p_ends.append(p_end)
            p = parts[-4]
            p_name = p[8:-2]
            p_names.append(p_name)
            sequence = str(record.seq)
            seqs.append(sequence)
            '''
        df = pd.DataFrame({'ID': ids, 'P_Name': p_names, 'P_start':p_starts, 'P_end':p_ends, 'Sequence': seqs})
        #print(df)
        print("############## Proteins for %s obtained ##############", x)
        #print(record.annotations["structured_comment"])
        handle.close()
        return df
    except Exception as e:
        print(e)

def get_feature_start_end(feature_location):
    start = None
    end = None

    for location in feature_location.parts:
        if location.start is None or location.end is None:
            continue
        if start is None:
            start = location.start
            end = location.end
        else:
            if location.start < start and location.strand == -1:
                start = location.start
            elif location.start > end and location.strand == 1:
                end = location.end
            if location.end > end and location.strand == -1:
                end = location.end
            elif location.end < start and location.strand == 1:
                start = location.start
    return start, end


def search_Entrez2(x):
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
                start, end = get_feature_start_end(feature.location)
                print(feature.location)
                print(start, end)
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
                    print(feature.qualifiers)
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
                print('!!!!! %s !!!!!' % feature.type)
                pl_ids.append(x)
                #### feature type
                p_type = feature.type
                p_types.append(p_type)
                #### location
                start, end = get_feature_start_end(feature.location)
                print(feature.location)
                print(start, end)
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
                    print(feature.qualifiers)
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
        print(df)
        print("############## Proteins for %s obtained ##############", x)
        # print(record.annotations["structured_comment"])
        handle.close()
    except Exception as e:
        print(e)


for plasmid in plasmids:
    search_Entrez2(plasmid)