#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 12:46:24 2024

@author: chris
"""

import os
import pandas as pd

os.chdir('/Users/chris/Desktop/prokka_dmnd/30qcov_30pid')


Results_top1 = pd.DataFrame(columns=['Genome', 'qseqid', 'sseqid', 'pident', 'length',
                                'mismatch', 'gapopen', 'qstart', 'qend', 'sstart',
                                'send', 'evalue', 'bitscore'])



genome_list = []

count=0
for file in os.listdir():
    if file.endswith('.csv'):
        
        blast_list = []
        genome = file.replace('non_mcy_hits_','')
        genome = genome.replace('.csv','')

        genome_list.append(genome)
        
        results = pd.read_csv(file, sep = ',', header=0)
        
        
        results['Genome'] = genome
        
        top_1 = results.groupby('gene').apply(lambda x: x.nlargest(1, 'bitscore')).reset_index(drop=True)
        top_1.to_csv('../all_top1/top1_'+genome+'.csv')
        Results_top1 = pd.concat([Results_top1, top_1])

Results_top1.to_csv('../all_top1/Compiled_top1_hits.csv', index=False)





