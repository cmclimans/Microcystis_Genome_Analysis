#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=12
#SBATCH --mem=2GB
#SBATCH --output=jobname_%J_stdout.txt
#SBATCH --error=jobname_%J_stderr.txt
#SBATCH --time=48:00:00
#SBATCH --job-name=muscle
#SBATCH --mail-user=cmclimans@ou.edu
#SBATCH --chdir=/scratch/cmclimans/Muscle
#
#################################################

source activate muscle

#muscle -align header_seqs_OnlyOrthos_processed_orthos_POACGBHM_03996.txt.csv.fa -clw test_out.aln

#./muscle5.1.linux_intel64 -in header_seqs_OnlyOrthos_processed_orthos_POACGBHM_03996.txt.csv.fa -clw -out aln_POACGBHM_03996.aln


#./muscle5.1.linux_intel64 -align header_seqs_OnlyOrthos_processed_orthos_POACGBHM_03996.txt.csv.fa -output aln.aln

while read line; do muscle -in $line -clwstrict -out FASconCAT/$line.aln; done < toalign.txt
