#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 12:38:18 2022

@author: chris
"""


# import os
import pandas as pd
#import sys
import argparse
parser = argparse.ArgumentParser()


parser.add_argument("-o", "--outfile", help="Outfile name for file with replaced headers")
parser.add_argument("-m", "--headers", help="Txt file with order of Headers to be repleaced")
parser.add_argument("-i", "--infile", help="Input sequence file (.fa)")
args = parser.parse_args()

fasta_file = args.infile
output_file_name = args.outfile



def main():
    # # Get genomes from a full GenID file if needed
    csv = pd.read_csv(args.headers, header = None)
    genomes = csv[0].tolist()
    
    
    # Read in a list of genomes from a single file
    Genome_order_list = []
    Genome_order_clean = []
    
    # Genome_order_list = pd.read_csv('rename_headers.csv', header = None).values.tolist()
    # for i in Genome_order_list:
    #     for j in i:
    #         if j not in Genome_order_clean:
    #             Genome_order_clean.append(j)
    
    
    if 'genomes' in locals():
        Genome_order_clean = genomes.copy()
        Genome_order_clean_tmp = Genome_order_clean.copy()
    elif Genome_order_clean in locals():
        Genome_order_clean_tmp = Genome_order_clean.copy()
    
    
    new_file = []
    Genome_order_clean_tmp = Genome_order_clean.copy()
    f = open(fasta_file, 'r')
    for line in f.readlines():        
        if line[0] == '>':    
            new_file.append(line.replace(line, '>'+Genome_order_clean_tmp[0].strip('.fa')+'\n'))
            Genome_order_clean_tmp.pop(0)
        else:
            new_file.append(line)
    f.close()
    newfileread = open(output_file_name,'w')
    for item in new_file:
        #newfileread.write(item+'\n')
        newfileread.write(item)
    newfileread.close()


main()
