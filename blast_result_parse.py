#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 07:49:32 2023

@author: Chris McLimans
@author_email: cmclimans@ou.edu

version = 0

This script is intended to parse blastn results (blast fmt 6) to identify gene presence in genomes.
This intends to check for complete gene sequences as well as full genes split around an assembly 
due to fragmented assembly. 

Currently this script is written to parse results for microcystin genes but will be developed
to search for other genes as long as gene legnths are provided and match qseqid. 

Goals:
    - Write section to identify complete (single full length >90%) matches of query [X]
    - Add section to confirm %ID of perfect match is acceptable [ ]
    - Write section to evaluate query with 2+ matches within single genome that start w/in 5% of qstart &| qend [ ]
        - Write section to check overlap and if split across multiple contigs [ ]
    - Write section to evaluate how much query aligns to genome for 2+ hits that do not start near qstart | qend [ ]
    - Convert to gene agnostic [ ]
    - Modify into functions [ ]
    - Convert to command line input [ ]
        - adjustable %ID and qlength acceptable options [ ]
    - Simplify repetitive sections and vectorize loops [ ]

"""


import os
import pandas as pd


os.chdir('/Users/chris/Desktop/all_mcy_blasts')
mcy_lengths = pd.read_csv('/Users/chris/Desktop/Karin_Botrys/mcy_blasts/mcy_lengths.csv', names = ['qseqid', 'Mcy_Length'])


Results = pd.DataFrame(columns=['Genome', 'qseqid', 'sseqid', 'pident', 'length',
                                'mismatch', 'gapopen', 'qstart', 'qend', 'sstart',
                                'send', 'evalue', 'bitscore'])

keep = ['Genome', 'qseqid', 'sseqid', 'pident', 'Q_aligned', 'qstart', 'qend',
        'Mcy_Length', 'sstart', 'send']

genome_list = []


for file in os.listdir('blast_all/'):
    if file.endswith('.txt'):
        genome = file.replace('out_','')
        genome = genome.replace('.txt','')
        genome_list.append(genome)
        results = pd.read_csv('blast_all/'+file, sep = '\t', names = ['qseqid', 'sseqid', 'pident',
                                                         'length', 'mismatch', 'gapopen',
                                                         'qstart', 'qend', 'sstart', 'send',
                                                         'evalue', 'bitscore'])
        
        
        results['Genome'] = genome
        results = pd.merge(results, mcy_lengths, on = 'qseqid')
        results['Mcy_Length'] = results['Mcy_Length'].astype(int)
        results['qend'] = results['qend'].astype(int)
        results['qstart'] = results['qstart'].astype(int)
        results['sstart'] = results['sstart'].astype(int)
        results['send'] = results['send'].astype(int)
        

        
        
        results['Q_aligned'] =  (results['qend'] - results['qstart'])
        results.loc[results['qstart'] == 1, 'Q_aligned'] += 1
        

                    
        
        results_keep = results[keep]
        results_keep = results_keep.sort_values(['Genome', 'qseqid', 'sstart'], ascending = True)  
        
        
        results_keep['Q_perc'] = results_keep['Q_aligned']/results_keep['Mcy_Length']*100
        
        results_small = results_keep[['Genome', 'qseqid', 'Q_perc', 'pident', 'qstart', 'qend', 'Mcy_Length', 'sseqid']]
        results_small.to_csv('blast_all_clean/out_'+genome+'.csv', index = False)

        
        Results = pd.concat([Results, results], ignore_index=True)


Genes = sorted(list(set(Results['qseqid'])))

# Start initial df to report gene presence
Result_df = pd.DataFrame(columns = Genes, index = genome_list)

# Count single matches that are at least 90% of query gene length as complete matches
## add in check for % ID
Perfect_matches = Results[abs(Results['Q_aligned'] - Results['Mcy_Length']) / Results['Mcy_Length'] <= 0.1]
Perfect_matches = Perfect_matches[['Genome', 'qseqid']]
Pmatch_indicies = Perfect_matches.index.to_list()

for index in Pmatch_indicies:
    tmpgene = Perfect_matches.loc[index]['qseqid']
    tmpgenome = Perfect_matches.loc[index]['Genome']
    Result_df.at[Result_df.index[Result_df[tmpgene].index == tmpgenome][0], tmpgene] = 'Yes'


Results['length'] = Results['length'].astype(int)
Results['mismatch'] = Results['mismatch'].astype(int)
Results['gapopen'] = Results['gapopen'].astype(int)
Results['qstart'] = Results['qstart'].astype(int)
Results['qend'] = Results['qend'].astype(int)
Results['sstart'] = Results['sstart'].astype(int)
Results['send'] = Results['send'].astype(int)


Result_df.to_csv('Perfect_Results.csv')


# Check for partial matches due to genome fragmentation

for genome in genome_list:
    genome = 'Ma_PCC9807'
    tmp_genome_df = Results[Results['Genome'] == genome]
    tmp_genome_genes = sorted(list(set(tmp_genome_df['qseqid'])))
    
    for gene in tmp_genome_genes:
        gene = 'mcyA'
        if pd.isna(Result_df.loc[genome][gene]) == True:
            if len(tmp_genome_df[tmp_genome_df['qseqid'] == gene]) == 1:
                print('here')
                # check length and % of qlength
                
            elif len(tmp_genome_df[tmp_genome_df['qseqid'] == gene]) > 1:
                tmp_genome_gene = tmp_genome_df[tmp_genome_df['qseqid'] == gene]
                
                mcy_length = int(mcy_lengths[mcy_lengths['qseqid'] == gene]['Mcy_Length'].iloc[0])
                max_start = int(mcy_length*0.05)
                min_end = mcy_length - int(mcy_length*0.05)
                
                tmp_genome_gene[(tmp_genome_gene['qstart'] >= 1 & tmp_genome_gene['qstart'] <= max_start)]
                
    

# Check for multiple gene hits in single genome
# if yes:                    
            # if >1 full (85%) hit:
                    # check if on different contigs
                    # report
        
        # check if at least 1 hit at qstart within 5% of gene length and at least 1 at within 5% of qend gene length
            # if 1 each
                # check if the two have overlapping region
                    # if yes = full hit
                    # report if spans different contigs
                # if no, determine missing length
                    # check if span multiple contigs - report
                    # full hit if >75% present
        
        # if multiple hits that don't start close to beginning or close to end
            # determine length of query found (non-overlap length)
            # check if on multiple contigs
            # report % length of query found
            # if at least 75% present but perhaps fragmented






    
"""
Old Code



# records = []
# for genome in Genomes:
#     genome_records = []
#     genome_df = Merged_keep[Merged_keep['Genome'] == genome]
#     for gene in list(set(genome_df['qseqid'])):
#         genome_records = []

#         gene_df = genome_df[genome_df['qseqid'] == gene]
#         n_hits = len(gene_df)
        
#         genome_records.append(genome)
#         genome_records.append(gene)
#         genome_records.append(n_hits)
        
#         records.append(genome_records)
        

# Max n hits is 6        

tmp_gene = genome_df[genome_df['qseqid'] == 'mcyA']




if len(tmp_gene) == 1:
    if tmp_gene['Q_aligned'] == tmp_gene['Mcy_Length']:
        print('Full gene')
        
    # elif something:
    #     partial gene mechanics

if len(tmp_gene) > 1:    
    check_perfect = tmp_gene['Q_aligned'] == tmp_gene['Mcy_Length']
    if check_perfect.any():
        perfect_index = check_perfect.index[check_perfect].values[0]
        perfect_match = pd.DataFrame(tmp_gene.loc[perfect_index]).transpose()
        
    elif (tmp_gene['qstart'] == 1).any() and (tmp_gene['qend'] == tmp_gene['Mcy_Length']).any():
        indicies = tmp_gene.index.to_list()
        check_start = tmp_gene['qstart'] == 1
        check_start_index = check_start.index[check_start].values[0]
        indicies.remove(check_start_index)
        
        # Check if the next records also starts with 1
            # if not check if overlaps with the end of the first 
            # if not check if starts immediately after the first record
            # else start is truly start

        # Check if any ranges overlap completely in non-first
            # Keep the largest range, remove the record for the smaller hit

        true_first = True
        for index, row in tmp_gene.iterrows():   
            print(index)
            print(row)
            
            
            index = 14
            if tmp_gene.loc[index]['qstart'] < tmp_gene.loc[check_start_index]['end']:
                true_first = False
            elif tmp_gene.loc[index]['qstart'] == tmp_gene.loc[check_start_index]['end']+1:
                print('second')
        
        if true_first:
            # Store first as the first record
                
        
    
#Merged_keep.to_csv('parsed_results.csv', index=False)      
"""

    