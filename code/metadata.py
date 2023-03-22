# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 20:52:42 2022

@author: sivan
"""

import pandas as pd
import logging
import re
import numpy as np
import sys
import getopt

#Read BLASTn file
plasmid_blast = pd.read_csv ("/path/to/Blast_ratio.csv", sep= ",",  header = 0)
plasmid_blast_copy = plasmid_blast.copy()

#Creating dataframe with unique plasmids ids
plasmid_blast_copy2 = plasmid_blast_copy[["qseqid"]]
plasmid_blast_copy2 = plasmid_blast_copy2.drop_duplicates()

#Read PLSDB file - chosen columns with taxonomy
with open("/path/to/the/plsdb/metadata/plsdb.tsv") as f:
    features_train = pd.read_csv(f, delimiter='\t', header=None,  error_bad_lines=False ,usecols=[1,28 ,30,32 , 34, 36, 38, 40])
    features_train_copy = features_train.copy()

#Merge by id the files in order to get a file with ids from the BLASTn results and taxonomy from PLSDB metadata    
plasmid_blast_copy2 = plasmid_blast_copy2.merge(features_train_copy, left_on= "0", right_on=1)

#Merge by id the files in order to get a file of BLASTn results and the taxonomy of thr plasmids 
plasmid_blast_copy = plasmid_blast_copy.merge(plasmid_blast_copy2, left_on= "qseqid", right_on="0")

plasmid_blast_copy2.to_csv("unique_taxonomy_metadata.csv")

