#!/bin/bash
#
#SBATCH --partition=largejobs
#SBATCH --ntasks=260
#SBATCH --nodes=13
#SBATCH --ntasks-per-node=20
#SBATCH --mem=60GB
#SBATCH --output=jobname_%J_stdout.txt
#SBATCH --error=jobname_%J_stderr.txt
#SBATCH --time=02:00:00
#SBATCH --job-name=mdmcleaner_lgjob
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cmclimans@ou.edu
#SBATCH --chdir=/scratch/cmclimans/pangenome/raw_fna
#
#################################################

#if anaconda 
eval "$(conda shell.bash hook)"

source activate mdmcleaner

#mdmcleaner makedb -o mdmcleaner_db_node/

mdmcleaner clean -i M60BS1.fna,88.fna -o ../decontaminated_fna_lgjob -t 200
