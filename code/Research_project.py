# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 13:30:50 2021

@author: sivan
"""

import Bio.Blast.Applications
import os,sys

#Turn spacer fasta file in to data base
from Bio.Blast.Applications import NcbimakeblastdbCommandline
makedb = NcbimakeblastdbCommandline(cmd='makeblastdb', dbtype='nucl', input_file=sys.argv[1], out = "spacersdb")
os.system(str(makedb))
stdout, stderr = makedb()

#Preform BLASTn on PLSDB
#This process was done on marine plasmids as well
from Bio.Blast.Applications import NcbiblastnCommandline
makeblast = NcbiblastnCommandline(query=sys.argv[2], db='/path/to/your/database/direcctory/spacersdb' ,evalue=0.00001, perc_identity= 90, ungapped=True, outfmt = "6 qseqid stitle sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore", out = "Blast.csv")
os.system(str(makeblast))
stdout, stderr = makeblast()



