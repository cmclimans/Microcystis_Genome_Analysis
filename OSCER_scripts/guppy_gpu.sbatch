#!/bin/bash
#
#SBATCH --partition=gpu_a100
#SBATCH --ntasks=10
#SBATCH --mem=10GB
#SBATCH --output=jobname_%J_stdout.txt
#SBATCH --error=jobname_%J_stderr.txt
#SBATCH --time=01:00:00
#SBATCH --job-name=guppy
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cmclimans@ou.edu
#SBATCH --chdir=/scratch/cmclimans/nanopore/Libs/20231107_1356_MN40502_FAT58221_5b44a2a9/
#
#################################################

#if anaconda 
eval "$(conda shell.bash hook)"


export PATH=$PATH:/home/cmclimans/software/ont-guppy/bin

module load CUDA/7.5.18
source activate guppy


#guppy_basecaller -i fast5_skip \
#-s fast5_gpu \
#-c dna_r9.4.1_450bps_sup_gpu.cfg -x 'auto'


guppy_barcoder -i fast5_gpu/pass/ -s barcodes_out --barcode_kits "SQK-RPB004" --num_barcoding_threads 10
