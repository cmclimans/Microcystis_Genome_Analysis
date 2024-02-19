#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 14:52:15 2023

@author: chris
"""

import sys
import pandas as pd
import os
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
parser = argparse.ArgumentParser()
import subprocess


parser.add_argument("-o", "--outfile", help="Outfile name for file with reverse complemented sequences as needed")
parser.add_argument("-p", "--primerID", help="acrB, glnA, helY, sbcC, amtB, trmR, lgt, dnaJ, ppnk, crp, murF, dinG, DUF6391, dICT, RmuC")
parser.add_argument("-i", "--indir", help="Input sequence directory")
parser.add_argument("-m", "--markerID", help="ID_file for marker IDs")
parser.add_argument("-t", "--threads", help="number of threads")


if len(sys.argv) < 3:
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()


seq_dict = {}
seq_dict.update({"acrB":Seq("TGGAG"),
                 "glnA":Seq("CATCG"),
                 "amtB":Seq("TCTTC"),
                 "sbcC":Seq("AACGC"),
                 "helY":Seq("GGCAA"),
                 "trmR":Seq("CGTTC"),
                 "lgt":Seq("GGACG"),
                 "dnaJ":Seq("GCAAA"),
                 "ppnk":Seq("GCCTC"),
                 "crp":Seq("GGTTC"),
                 "murF":Seq("CCTGC"),
                 "dinG":Seq("CGGGA"),
                 "DUF6391":Seq("GCAGT"),
                 "dICT":Seq("AGTCT"),
                 "RmuC":Seq("GATTT")
                 })



GCA_key = pd.read_csv(args.markerID, names = ['ID', 'Abrev'])
GCA_dict = GCA_key.set_index('ID').to_dict()['Abrev']


startseq = seq_dict[args.primerID]
records = []
missing_seqs = []
has_seqs = []
count = 0

for file in os.listdir(args.indir):
    if file.endswith('.fna'):
                
        # GCA ID
        if 'GC' in file:
            ID = file.split('_GC')[1]
            ID = 'GC'+ID
            ID = ID.split('.fna')[0]
        
        # Botrys ID
        elif 'S5' in file:
            ID = file.split('_S5')[1]
            ID = 'S5'+ID
            ID = ID.split('.fna')[0]
            
        if ID in GCA_dict:
            ID = GCA_dict[ID]


        #record = SeqIO.parse(open(args.indir+'/'+file),"fasta")
        
        
        for record in SeqIO.parse(open(args.indir+'/'+file),"fasta"):
            if len(record.seq) > 0:
                if record.seq[0:5] == startseq:
                    newrecord = record.reverse_complement()
                    newrecord.id = ID
                    newrecord.description = ''
                    records.append(newrecord)
                else:
                    newrecord = record
                    newrecord.id = ID
                    newrecord.description = ''
                    records.append(newrecord)
                
                has_seqs.append(ID)
                
            else:
                missing_seqs.append(ID)
                
        


SeqIO.write(records, f'{args.indir}/{args.outfile}.fasta', "fasta")

if len(missing_seqs) > 0:
    with open(f'{args.indir}/{args.primerID}_missing_seqs.txt','w') as f:
              f.writelines([f"{item}\n" for item in missing_seqs])

if len(has_seqs) > 0:
    with open(f'{args.indir}/{args.primerID}_seq_headers.txt','w') as f:
              f.writelines([f"{item}\n" for item in has_seqs])






# # Dev area for seq processing

# file = 'test.fna'
# startseq = 'CCTTC's
# records = []
# for record in SeqIO.parse(open('/Users/chris/Desktop/test.fna'),"fasta"):
#     if len(record.seq) > 0:
#         if record.seq[0:5] == startseq:
#             newrecord = record.reverse_complement()
#             newrecord.id = ID
#             newrecord.description = ''
#             records.append(newrecord)
#     else:
#         print(file.strip('.fna'))




# primer = 'acrB'
# infile = '/Users/chris/Desktop/acrB_combo.fna'
# outfile = f'/Users/chris/Desktop/{primer}_align.aln'

muscle_cmd = f'source activate muscle; muscle -in {args.indir}/{args.outfile}.fasta -fasta -out {args.indir}/{args.primerID}_align.aln'
alignment = subprocess.run(muscle_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

if alignment.returncode != 0:
    print(f'Error occurred: {alignment.stderr}')
    exit(1)


if os.path.isdir(f'trees/{args.primerID}') is False:
    os.mkdir(f'trees/{args.primerID}')




tree_cmd1 = f'module load RAxML/8.2.11-intel-2018a-pthreads-avx2; raxmlHPC -T {args.threads} \
    -o Ma_SC_T_19800800_S464 -s {args.indir}/{args.primerID}_align.aln \
    -n {args.primerID} -w /scratch/cmclimans/sequences/trees/{args.primerID} \
    -m GTRCAT -b 522 -p 22 -N 100; \
    raxmlHPC -T {args.threads} -m GTRCAT -p 22 -# 100 -s {args.indir}/{args.primerID}_align.aln \
    -n {args.primerID}_ML -w /scratch/cmclimans/sequences/trees/{args.primerID} -o Ma_SC_T_19800800_S464 \
    raxmlHPC -T {args.threads} -m GTRCAT -p 22 -f b -t trees/{args.primerID}/RAxML_bestTree.{args.primerID}_ML; \
    -z trees/{args.primerID}/RAxML_bootstrap.{args.primerID} -n {args.primerID}_bootstrap_besttree_out -o Ma_SC_T_19800800_S464 \
    -w /scratch/cmclimans/sequences/trees/{args.primerID}'
    
subprocess.run(tree_cmd1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


if os.path.isdir(f'/scratch/cmclimans/sequences/trees/{args.primerID}/tmp_files') is False:
    os.mkdir(f'/scratch/cmclimans/sequences/trees/{args.primerID}/tmp_files')

os.chdir(f'/scratch/cmclimans/sequences/trees/{args.primerID}')

for file in os.listdir():
    if 'RUN' in file:
        shutil.move(file, 'tmp_files')
        
        
if os.path.isdir('/scratch/cmclimans/sequences/best_trees') is False:
    os.mkdir('/scratch/cmclimans/sequences/best_trees')


for file in os.listdir():
    if 'bootstrap_besttree_out' in file:
        shutil.copy(file, '../../best_trees')
       



    


# # Dev area for RAxML
# tree_cmd2 = f'raxmlHPC -T {args.threads} -m GTRCAT -p 22 -# 100 -s {args.indir}/{args.primerID}_align.aln \
#     -n {args.primerID}_ML -w /scratch/cmclimans/sequences/trees/{args.primerID} -o Ma_SC_T_19800800_S464 \
#         raxmlHPC -T {args.threads} -m GTRCAT -p 22 -f b -t trees/{args.primerID}/RAxML_bestTree.{args.primerID}_ML'
    
# tree_cmd3 = f'raxmlHPC -T {args.threads} -m GTRCAT -p 22 -f b -t trees/{args.primerID}/RAxML_bestTree.{args.primerID}_ML \
#     -z trees/{args.primerID}/RAxML_bootstrap.{args.primerID} -n {args.primerID}_bootstrap_besttree_out -o Ma_SC_T_19800800_S464 \
#         -w /scratch/cmclimans/sequences/trees/{args.primerID}'


#subprocess.run(tree_cmd2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
#subprocess.run(tree_cmd3, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


# if tree.returncode != 0:
#     print(f'Error occurred: {alignment.stderr}')
#     print(f'Error occurred: {alignment.stdout}')

#     exit(1)

    
    
    
    
    
