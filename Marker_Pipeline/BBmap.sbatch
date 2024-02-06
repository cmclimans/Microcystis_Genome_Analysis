#!/bin/bash
#
#SBATCH --partition=debug
#SBATCH --ntasks=4
#SBATCH --mem=16GB
#SBATCH --output=jobname_%J_stdout.txt
#SBATCH --error=jobname_%J_stderr.txt
#SBATCH --time=00:30:00
#SBATCH --job-name=BBmapF
#SBATCH --mail-user=cmclimans@ou.edu
#SBATCH --chdir=/scratch/cmclimans/botrys
#
#################################################

module load BBMap

# glnA_F
#while read line; do msa.sh -Xmx16g in=fna/$line out=BBmap_out_glnA/BBmap_glnA_F_$line literal=ATGCCCGAAACGCCACAAGAAGTC rcomp=t cutoff=0.875; done < fna/files.txt

# helY_F
#while read line; do msa.sh -Xmx16g in=fna/$line out=BBmap_out_helY/BBmap_helY_F_$line literal=GCTCTTTCTAACCAAAAATTCCGGG rcomp=t cutoff=0.875; done < fna/files.txt

# acrB_F
#while read line; do msa.sh -Xmx16g in=fna/$line out=BBmap_out_acrB/BBmap_acrB_F_$line literal=CCTTCTTTTAGTATCGATCAGAGAG rcomp=t cutoff=0.875; done < fna/files.txt


# amtB_F
#while read line; do msa.sh -Xmx16g in=fna/$line out=BBmap_out_amt/BBmap_amtB_F_$line literal=AACGAATTGCTCGAAAATCCGGGG rcomp=t cutoff=0.875; done < fna/files.txt


# sbcC_F
while read line; do msa.sh -Xmx16g in=fna/$line out=BBmap_out_sbcC/BBmap_sbcC_F_$line literal=CAAGCGGACGARTTTATGAAACAG rcomp=t cutoff=0.875; done < fna/files.txt
