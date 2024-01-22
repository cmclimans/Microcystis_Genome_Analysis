#!/bin/bash
#
#SBATCH --partition=debug
#SBATCH --ntasks=1
#SBATCH --mem=1GB
#SBATCH --output=jobname_%J_stdout.txt
#SBATCH --error=jobname_%J_stderr.txt
#SBATCH --time=00:30:00
#SBATCH --job-name=seqtk
#SBATCH --mail-user=cmclimans@ou.edu
#SBATCH --chdir=/home/cmclimans/
#
#################################################


source activate seqtk

while read line; do seqtk subseq /scratch/cmclimans/prep_analysis/peppan_out/peppan_out_preps.allele.fna /scratch/cmclimans/prep_analysis/single_geno_genes/reformatted/$line > /scratch/cmclimans/prep_analysis/single_geno_genes/reformatted/seqs_$line; done < /scratch/cmclimans/prep_analysis/single_geno_genes/reformatted/seq_list.txt


