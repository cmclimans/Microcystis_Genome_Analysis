#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 14:52:15 2023

@author: chris
"""

# import os
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
parser = argparse.ArgumentParser()


parser.add_argument("-o", "--outfile", help="Outfile name for file with reverse complemented sequences as needed")
parser.add_argument("-p", "--primerID", help="acrB, glnA, helY, sbcC, amtB")
parser.add_argument("-i", "--infile", help="Input sequence file (.fa)")
args = parser.parse_args()


seq_dict = {}
seq_dict.update({"acrB":Seq("TGGAGAT"),
                 "glnA":Seq("CATCGGG"),
                 "amtB":Seq("TCTTCGG"),
                 "sbcC":Seq("AACGCCA"),
                 "helY":Seq("GGCAATT")})




startseq = seq_dict[args.primerID]

records = []
count = 0
for record in SeqIO.parse(args.infile,"fasta"):
    ID = record.id
    type(record.id)
    # print(record.id)
    # print(record.seq)
    if record.seq[0:7] == startseq:
        newrecord = record.reverse_complement()
        newrecord.id = ID
        records.append(newrecord)
    else:
        records.append(record)

        
SeqIO.write(records, args.outfile, "fasta")


