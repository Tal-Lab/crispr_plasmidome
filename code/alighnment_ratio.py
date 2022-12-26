# -*- coding: utf-8 -*-
"""
Created on Sat May  7 22:01:30 2022

@author: sivan
"""


import pandas as pd
import sys

#Read files using sis.argv in order to summon the files names
plasmid_blast = pd.read_csv (sys.argv[1], sep= ",",  header = 0)
plasmid_blast_copy = plasmid_blast.copy()

#Finding ratio of spacer sequence length and alignment length
plasmid_blast_copy["ratio"] = plasmid_blast_copy ["length"]/plasmid_blast_copy ["spacer_length"]

#Create CSV using sis.argv in order to summon the output files names
plasmid_blast_copy.to_csv(sys.argv[2])