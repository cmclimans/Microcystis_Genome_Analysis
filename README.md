<img src="Bloom_Image.jpeg" width="700" height="500">

# Microcystis_Genome_Analysis  


## Project Description

This repository contains scripts and workflows required to perform genome analyses on _Microcystis_. 

## Project Goals
1. Collate and QC _Microcystis_ genomes  
	1a. Collate Genomes  
	- Place fasta files into a common directory.
		- Create a rename .csv file to rename fasta files to a common format.  
	- Confirm fasta files are all (.fna) extension to ensure downstream commands work. Rename with .fna if needed.
	- Rename fasta files, as needed, to a common format or scheme using rename_fasta.py script and rename .csv file.  

	1b. Decontaminate with mdmclener  
	- Install mdmcleaner via anaconda (https://github.com/KIT-IBG-5/mdmcleaner) or using  mdmcleaner.yml file  
		- If not installed, run `mdmcleaner makedb` to install databases  
		- specify output db directory with the -o option  
	- This will take a while (~14 hours on a 100 mbps network, per developer)  
	- run `mdmcleaner clean -i FASTA_FILES -o OUTPUT_DIR -t THREADS --overview_files_basename RUN_NAME` --no_filterfasta True`  
		- one output file per genome will be produced, we need the following files:  
			- <genome-name>_keep.fasta which contains safe, uncontaminated contigs  
			- <genome-name>_evaluate_low.fasta contains contigs which _may_ have contaminated reads  
			- <genome-name>_evaluate_high.fasta contains contigs likely to be contaminated  
			- <genome-name>_delete.fasta poor quality contigs likely with definete contamination which should be removed  

	1c. QC Decontamination with QUAST  

3. Pangenome analysis of genomes
4. Phylogenetic inference
