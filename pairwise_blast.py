#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 08:16:48 2024

@author: Chris McLimans

Purpose:
    - Following extraction of gene sequence from genomes and correcting for
      reverse complements and start position errors as needed, this script
      will run blast for all pairwise comparisons of sequenes to record
      pairwise %ID of sequences. 

Goals:
    - Read in headers from fasta to later iterate over [X]
    - Use each seq as query, then blast against the remaining sequences as 
      subject, report best hit per subject [X]
      - Remove query seq from 
    - add flexibility to switch nucl and aa []
    - For each query-subject blast, store results in table (use outfmt 6) [X]
    - (Maybe) reformat results table as pairwise matrix []
    - Report average %ID, worst %ID, best %ID, etc. as an output file [X]
    - Integrate to prior blast parse scripts outputs []

"""


import os
import pandas as pd
import subprocess
from Bio import SeqIO
#import shutil


if os.path.isdir('aln_tmp/') is False:
   os.mkdir('aln_tmp/')
    

# Import the gene sequences for subsequent use

Blast_results = pd.DataFrame()
allStats = pd.DataFrame()

genes = ['mcyA', 'mcyB', 'mcyC', 'mcyD', 'mcyE', 'mcyF', 'mcyG', 'mcyH', 'mcyI', 'mcyJ']

for gene in genes:
    
    gene = 'mcyA'
    Blast_results = pd.DataFrame()
    os.chdir('/Users/chris/Desktop/alignments/'+gene+'_align')
    all_seqs = SeqIO.index(gene+'_fixed.faa', 'fasta')
    
    # Store the sequence header names to iterate over
    seqs = []
    for record in SeqIO.parse(gene+'_fixed.faa', 'fasta'):
        seqs.append(record.id)
    seqs_copy = seqs.copy()
    
    
    # Extract the sequence to use as the query sequence
    count = 0
    no_hits = []
    for sequence in seqs:
        if sequence != seqs[len(seqs)-1]:
            query = all_seqs.get_raw(sequence).decode()
            count +=1
            
            print(gene+', '+sequence+', '+str(count)+'/'+str(len(seqs)))
            
            with open ('aln_tmp/tmp_query.faa', 'w') as f:
                f.write(query)
            
            # Remove the query seq ID from the list of sequences
            seqs_copy.remove(sequence)
            
            # Store fasta sequences that are not the query sequence
            records = [] 
            for sequence_sub in seqs_copy:
                records.append(all_seqs.get_raw(sequence_sub).decode())
            
            # Convert the non-query sequences to subject file
            record_len = len(records)
            subjects = '\n'.join([''.join(records[0:record_len])])
            with open('aln_tmp/tmp_subjects.faa', 'w') as f:
                f.write(subjects)
            
            
            query = '/Users/chris/Desktop/alignments/'+gene+'_align/aln_tmp/tmp_query.faa'
            subjets = '/Users/chris/Desktop/alignments/'+gene+'_align/aln_tmp/tmp_subjects.faa'
            
            blast_cmd = f'source activate blast; blastp -query {query} -subject {subjets} -outfmt 6 -subject_besthit'
            result = subprocess.run(blast_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            
            if result.returncode != 0:
                print(f'Error occurred: {result.stderr}')
                exit(1)
            else:
                blast_output = result.stdout
                
            rows = blast_output.split('\n')
            table = [row.split('\t') for row in rows]
            if len(table) > 1:
                blast_table = pd.DataFrame(table,columns=['qseqid', 'sseqid', 'pident', 'length',
                                                'mismatch', 'gapopen', 'qstart', 'qend', 'sstart',
                                                'send', 'evalue', 'bitscore'])
                
            elif len(table) <= 1:
                no_hits.append(sequence)
            
            Blast_results = pd.concat([Blast_results,blast_table])
        
        # if os.path.isdir('aln_tmp/') is True:
        #     shutil.rmtree('aln_tmp')
            
        Blast_results.columns =['qseqid', 'sseqid', 'pident', 'length',
                                        'mismatch', 'gapopen', 'qstart', 'qend', 'sstart',
                                        'send', 'evalue', 'bitscore']
        
        Blast_results.to_csv(gene+'_pairwise_prot.csv', index=False)
        
        Blast_results['pident'] = Blast_results['pident'].astype(float)
        
        
        with open('no_hits.txt', 'w') as f:
            for item in no_hits:
                f.write(item)
        
        
        
        Blast_results = pd.read_csv('/Users/chris/Desktop/alignments/'+gene+'_align/'+gene+'_pairwise_prot.csv')
        
    
        stats = {}
        stats['Max'] = Blast_results['pident'].max()
        stats['Min'] = Blast_results['pident'].min()
        stats['Mean'] = Blast_results['pident'].mean()
        stats['StDev'] = Blast_results['pident'].std()
        stats['CompleteSequences'] = len(list(set(Blast_results['sseqid'])))
        
        Stats_DF = pd.DataFrame.from_dict(stats, orient='index', columns = [gene]).T
        Stats_DF.to_csv('/Users/chris/Desktop/alignments/'+gene+'_align/Stats_prot.csv', index=False)
    
        allStats = pd.concat([allStats, Stats_DF])
    
allStats.to_csv('/Users/chris/Desktop/alignments/blast_stats_prot.csv')



