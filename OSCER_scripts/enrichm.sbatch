#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=6
#SBATCH --mem=8GB
#SBATCH --output=jobname_%J_stdout.txt
#SBATCH --error=jobname_%J_stderr.txt
#SBATCH --time=48:00:00
#SBATCH --job-name=enrichm
#SBATCH --mail-user=cmclimans@ou.edu
#SBATCH --chdir=/scratch/cmclimans/
#
#################################################

source activate enrichm
export ENRICHM_DB=/scratch/cmclimans/enrichm_db

#enrichm annotate --ko --output all_genomes_enrichmAnnotated/ --genome_directory all_genomes/ --threads 12 --force



#enrichm annotate --ko --output all_genomes_enrichmAnnotated --genome_files GCA_000010625.1_ASM1062v1_genomic.fna --threads 8 --force

#enrichm classify --genome_and_annotation_matrix all_genomes_enrichmAnnotated/ko_frequency_table.tsv --output all_genomes_enrichmAnnotated_classify --aggregate --force

enrichm enrichment --annotate_output all_genomes_enrichmAnnotated/ --metadata enrichment_groups_troph.csv --ko --processes 6 --output all_genomes_enrichmAnnotated_enrichment --force

#enrichm enrichment --annotate_output all_genomes_enrichmAnnotated/ --output all_genomes_enrichmAnnotated_enrichment_debug/ --metadata enrichment_groups.csv --ko --processes 6 --force

