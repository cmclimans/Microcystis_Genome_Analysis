#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=12
#SBATCH --mem=24GB
#SBATCH --output=jobname_%J_stdout.txt
#SBATCH --error=jobname_%J_stderr.txt
#SBATCH --time=12:00:00
#SBATCH --job-name=checkm
#SBATCH --mail-user=cmclimans@ou.edu
#SBATCH --chdir=/home/cmclimans/122_Genome_FNA
#
#################################################
#hostname

module load CheckM

mkdir leftover_checkm_out
checkm lineage_wf leftover/ leftover_checkm_out/ -t 12


