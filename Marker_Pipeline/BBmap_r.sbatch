#!/bin/bash
#
#SBATCH --partition=debug
#SBATCH --ntasks=4
#SBATCH --mem=16GB
#SBATCH --output=jobname_%J_stdout.txt
#SBATCH --error=jobname_%J_stderr.txt
#SBATCH --time=00:30:00
#SBATCH --job-name=BBmapR
#SBATCH --mail-user=cmclimans@ou.edu
#SBATCH --chdir=/scratch/cmclimans/botrys
#
#################################################

module load BBMap

# glnA_R
#while read line; do msa.sh -Xmx16g in=fna/$line out=BBmap_out_glnA/BBmap_glnA_R_$line literal=CATCGGGTTAACTTCGTTATCGAG rcomp=t cutoff=0.875; done < fna/files.txt

# helY_R
#while read line; do msa.sh -Xmx16g in=fna/$line out=BBmap_out_helY/BBmap_helY_R_$line literal=GGCAATTTAGATGTATTTCCTGTC rcomp=t cutoff=0.875; done < fna/files.txt

# acrB_R
#while read line; do msa.sh -Xmx16g in=fna/$line out=BBmap_out_acrB/BBmap_acrB_R_$line literal=TGGAGATTATCGGAGGCAATATA rcomp=t cutoff=0.875; done < fna/files.txt


# amtB_R
#while read line; do msa.sh -Xmx16g in=fna/$line out=BBmap_out_amt/BBmap_amtB_R_$line literal=TCTTCGGGAGTAACACGAATACCG rcomp=t cutoff=0.875; done < fna/files.txt


# sbcC_R
while read line; do msa.sh -Xmx16g in=fna/$line out=BBmap_out_sbcC/BBmap_sbcC_R_$line literal=AACGCCATCTTTTTTCTTCTCCTC rcomp=t cutoff=0.875; done < fna/files.txt

