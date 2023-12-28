#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 12:38:18 2022

@author: chris
"""


import pandas as pd
#import sys
import argparse
parser = argparse.ArgumentParser()


parser.add_argument("-o", "--outfile", help="Outfile name for file with samples removed that don't have a primer hit")
parser.add_argument("-p1", "--primer1", help="first primer sequences check file")
parser.add_argument("-p2", "--primer2", help="second primer sequences check file")
parser.add_argument("-i", "--infile", help="Input text file from list of starting fasta list")
args = parser.parse_args()



    

def split(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])


#dd
def main():
    # # Example using args
    sample_to_remove = []
    fasta_list = pd.read_csv(args.infile, header = None)  
    
    p1_list = pd.read_csv(args.primer1, sep=":", header= None)
    p2_list = pd.read_csv(args.primer2, sep=":", header= None)
    
    
    p1_tmp_ID = p1_list[0][1]
    p1_str = split(p1_tmp_ID, '_', 2)[0]+'_' 
    p2_tmp_ID = p2_list[0][1]
    p2_str = split(p2_tmp_ID, '_', 2)[0]+'_'
    
    
    p1_missing = p1_list.loc[p1_list[1] == 0, [0]].values.tolist()
    p2_missing = p2_list.loc[p2_list[1] == 0, [0]].values.tolist()

    
    for item in p1_missing:
        for ID in item:
            sample_to_remove.append(ID.replace(p1_str,''))
    for item in p2_missing:
        for ID in item:
            sample_to_remove.append(ID.replace(p2_str,''))

    
    fasta_list_updated = fasta_list[~fasta_list[0].isin(sample_to_remove)]
    fasta_list_updated.to_csv(args.outfile, index = False, header=None)

main()
