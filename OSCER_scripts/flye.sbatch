#!/bin/bash
#
#SBATCH --partition=debug
#SBATCH --ntasks=20
#SBATCH --mem=10GB
#SBATCH --output=jobname_%J_stdout.txt
#SBATCH --error=jobname_%J_stderr.txt
#SBATCH --time=00:30:00
#SBATCH --job-name=flye
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cmclimans@ou.edu
#SBATCH --chdir=/scratch/cmclimans/nanopore
#
#################################################

#if anaconda 
eval "$(conda shell.bash hook)"


source activate flye

#mkdir B1_assembly
#flye --nano-hq B1.fastq --out-dir B1_assembly --threads 20 --iterations 3 --scaffold


#mkdir B2_assembly
#flye --nano-hq B2.fastq	--out-dir B2_assembly --threads	20 --iterations	3 --scaffold


mkdir B6_assembly_meta
flye --nano-hq B6.fastq	--out-dir B6_assembly_meta --threads 20 --iterations 3 --scaffold --meta


#mkdir B7_assembly
#flye --nano-hq B7.fastq --out-dir B7_assembly --threads	20 --iterations	3 --scaffold


#mkdir B10_assembly
#flye --nano-hq B10.fastq --out-dir B10_assembly --threads 20 --iterations 3 --scaffold


 
