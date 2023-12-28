#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=12
#SBATCH --mem=20GB
#SBATCH --output=jobname_%J_stdout.txt
#SBATCH --error=jobname_%J_stderr.txt
#SBATCH --time=12:00:00
#SBATCH --job-name=PrimerTree
#SBATCH --mail-user=cmclimans@ou.edu
#SBATCH --chdir=/scratch/cmclimans/botrys
#
#################################################


# Logistics to set working marker directory
if [[ -z "$1" ]]; then
    echo "Must provide working directory name in command line position 1" 1>&2
    exit 1

if [[ "$1" ]]; then
    echo "Will perform analysis and write files to $1" 1>&2
fi

if [ -d "$1" ]
then
    echo "Directory /path/to/dir exists."
    exit 1
fi



# Create a new Directory for the output of the analysis
# positonal argument from command line. Place the directory ID as the first position argument
mkdir $1


# Put the fasta format genomes into a directory called fna
# Make sure the list of files is updated and ready for analysis
ls fna/*.fa* > fna/files.txt


# Load BBMap software (alter as needed for installation)
# this works on OSCER
module load BBMap


### Step 1
# Primer Coordinate Extraction
# Use BBMap to find the position in each genome that the forward and reverse primers aligns.
# This outputs a SAM file that can then be used by BBMAp to extract the sequences including/between primers.
# need -Xmx20g to specify 20g memory to Java (else large genomes fail), rcomp needed becuase genomes are not always in the same orientation
# cutoff is a percent here allowing three mismatches given quality of MAGs, however hits are generally 100% but allowing flexibility
#
#
###

while read line; do msa.sh -Xmx20g in=fna/$line out=$1/BBmap_out_glnA/BBmap_glnA_F_$line literal=ATGCCCGAAACGCCACAAGAAGTC rcomp=t cutoff=0.875; done < fna/files.txt
while read line; do msa.sh -Xmx20g in=fna/$line out=$1/BBmap_out_glnA/BBmap_glnA_R_$line literal=CATCGGGTTAACTTCGTTATCGAG rcomp=t cutoff=0.875; done < fna/files.txt

while read line; do msa.sh -Xmx20g in=fna/$line out=$1/BBmap_out_helY/BBmap_helY_F_$line literal=GCTCTTTCTAACCAAAAATTCCGGG rcomp=t cutoff=0.875; done < fna/files.txt
while read line; do msa.sh -Xmx20g in=fna/$line out=$1/BBmap_out_helY/BBmap_helY_R_$line literal=GGCAATTTAGATGTATTTCCTGTC rcomp=t cutoff=0.875; done < fna/files.txt

while read line; do msa.sh -Xmx20g in=fna/$line out=$1/BBmap_out_acrB/BBmap_acrB_F_$line literal=CCTTCTTTTAGTATCGATCAGAGAG rcomp=t cutoff=0.875; done < fna/files.txt
while read line; do msa.sh -Xmx20g in=fna/$line out=$1/BBmap_out_acrB/BBmap_acrB_R_$line literal=TGGAGATTATCGGAGGCAATATA rcomp=t cutoff=0.875; done < fna/files.txt

while read line; do msa.sh -Xmx20g in=fna/$line out=$1/BBmap_out_amt/BBmap_amtB_F_$line literal=AACGAATTGCTCGAAAATCCGGGG rcomp=t cutoff=0.875; done < fna/files.txt
while read line; do msa.sh -Xmx20g in=fna/$line out=$1/BBmap_out_amt/BBmap_amtB_R_$line literal=TCTTCGGGAGTAACACGAATACCG rcomp=t cutoff=0.875; done < fna/files.txt

while read line; do msa.sh -Xmx20g in=fna/$line out=$1/BBmap_out_sbcC/BBmap_sbcC_F_$line literal=CAAGCGGACGARTTTATGAAACAG rcomp=t cutoff=0.875; done < fna/files.txt
while read line; do msa.sh -Xmx20g in=fna/$line out=$1/BBmap_out_sbcC/BBmap_sbcC_R_$line literal=AACGCCATCTTTTTTCTTCTCCTC rcomp=t cutoff=0.875; done < fna/files.txt



### Step 2
# Extract Primer Seqs
# Use the coordinates of the primers from each genome as identified in Step 1 to 'cut' the sequences between
# the primers (inclues the primer sequence)
# sam1 and sam2 specify the F and R SAM files for a sample
# include=t includes the primer sequence in the extracted sequence, else only get the sequence between the primers
# fake=f ensures that empty files are written IF the primer does not align, else writes a fasta with a header from the source genome
#
# Then the the Reverse_sequences_forPrimer.py python scrpt written by Chris McLimans to check the orientation of each sequence
# since some genomes may be 5 -> 3 and some 3 -> 5 from assembly, the extracted sequences may be different orientations.
# this script will correct them so they're all 5 -> 3 and will align properly later.
#
###

# Move into the Marker working directory
cd $1

# Make directories for each primer pairing
mkdir BBmap_out_acrB_glnA
mkdir BBmap_out_amtB_sbcC
mkdir BBmap_out_glnA_helY
mkdir BBmap_out_helY_sbcC

# Copy the file name list of the fasta files to the primer pairing directories (needed later to
# rename the extracted sequences)
cp ../fna/files.txt BBmap_out_acrB_glnA/
cp ../fna/files.txt BBmap_out_amtB_sbcC/
cp ../fna/files.txt BBmap_out_glnA_helY/
cp ../fna/files.txt BBmap_out_helY_sbcC/


# Extract the sequences between and including the primer for each primer pairing
# Writes ONE fasta per sample regardless of whether a primer hits or not.
while read line; do cutprimers.sh in=../fna/$line out=BBmap_out_acrB_glnA/seq_acrB_$line sam1=BBmap_out_acrB/BBmap_acrB_F_$line sam2=BBmap_out_acrB/BBmap_acrB_R_$line include=t fake=f; done < BBmap_out_acrB_glnA/files.txt
while read line; do cutprimers.sh in=../fna/$line out=BBmap_out_acrB_glnA/seq_glnA_$line sam1=BBmap_out_glnA/BBmap_glnA_F_$line sam2=BBmap_out_glnA/BBmap_glnA_R_$line include=t fake=f; done < BBmap_out_acrB_glnA/files.txt

while read line; do python Reverse_sequences_forPrimer.py -i BBmap_out_acrB_glnA/seq_acrB_$line -p 'acrB' -o BBmap_out_acrB_glnA/oriented_seq_acrB_$line; done < BBmap_out_acrB_glnA/files.txt
while read line; do python Reverse_sequences_forPrimer.py -i BBmap_out_acrB_glnA/seq_glnA_$line -p 'glnA' -o BBmap_out_acrB_glnA/oriented_seq_glnA_$line; done < BBmap_out_acrB_glnA/files.txt



while read line; do cutprimers.sh in=../fna/$line out=BBmap_out_amtB_sbcC/seq_amtB_$line sam1=BBmap_out_amt/BBmap_amtB_F_$line sam2=BBmap_out_amt/BBmap_amtB_R_$line include=t fake=f; done < BBmap_out_amtB_sbcC/files.txt
while read line; do cutprimers.sh in=../fna/$line out=BBmap_out_amtB_sbcC/seq_sbcC_$line sam1=BBmap_out_sbcC/BBmap_sbcC_F_$line sam2=BBmap_out_sbcC/BBmap_sbcC_R_$line include=t fake=f; done < BBmap_out_amtB_sbcC/files.txt

while read line; do python Reverse_sequences_forPrimer.py -i BBmap_out_amtB_sbcC/seq_amtB_$line -p 'amtB' -o BBmap_out_amtB_sbcC/oriented_seq_amtB_$line; done < BBmap_out_amtB_sbcC/files.txt
while read line; do python Reverse_sequences_forPrimer.py -i BBmap_out_amtB_sbcC/seq_sbcC_$line -p 'sbcC' -o BBmap_out_amtB_sbcC/oriented_seq_sbcC_$line; done < BBmap_out_amtB_sbcC/files.txt



while read line; do cutprimers.sh in=../fna/$line out=BBmap_out_glnA_helY/seq_glnA_$line sam1=BBmap_out_glnA/BBmap_glnA_F_$line sam2=BBmap_out_glnA/BBmap_glnA_R_$line include=t fake=f; done < BBmap_out_glnA_helY/files.txt
while read line; do cutprimers.sh in=../fna/$line out=BBmap_out_glnA_helY/seq_helY_$line sam1=BBmap_out_helY/BBmap_helY_F_$line sam2=BBmap_out_helY/BBmap_helY_R_$line include=t fake=f; done < BBmap_out_glnA_helY/files.txt

while read line; do python Reverse_sequences_forPrimer.py -i BBmap_out_glnA_helY/seq_glnA_$line -p 'amtB' -o BBmap_out_glnA_helY/oriented_seq_glnA_$line; done < BBmap_out_glnA_helY/files.txt
while read line; do python Reverse_sequences_forPrimer.py -i BBmap_out_glnA_helY/seq_helY_$line -p 'sbcC' -o BBmap_out_glnA_helY/oriented_seq_helY_$line; done < BBmap_out_glnA_helY/files.txt



while read line; do cutprimers.sh in=../fna/$line out=BBmap_out_helY_sbcC/seq_helY_$line sam1=BBmap_out_helY/BBmap_helY_F_$line sam2=BBmap_out_helY/BBmap_helY_R_$line include=t fake=f; done < BBmap_out_helY_sbcC/files.txt
while read line; do cutprimers.sh in=../fna/$line out=BBmap_out_helY_sbcC/seq_sbcC_$line sam1=BBmap_out_sbcC/BBmap_sbcC_F_$line sam2=BBmap_out_sbcC/BBmap_sbcC_R_$line include=t fake=f; done < BBmap_out_helY_sbcC/files.txt

while read line; do python Reverse_sequences_forPrimer.py -i BBmap_out_helY_sbcC/seq_helY_$line -p 'helY' -o BBmap_out_helY_sbcC/oriented_seq_helY_$line; done < BBmap_out_helY_sbcC/files.txt
while read line; do python Reverse_sequences_forPrimer.py -i BBmap_out_helY_sbcC/seq_sbcC_$line -p 'sbcC' -o BBmap_out_helY_sbcC/oriented_seq_sbcC_$line; done < BBmap_out_helY_sbcC/files.txt




### Step 3
# Check Exracted Primers
# uses grep to check which each sample has an extracted sequence for each primer pair
# writes text file in the format of SAMPLE:1 where 1 is the number of sequences extracted for the SAMPLE
###



cd BBmap_out_acrB_glnA
grep -c '>' oriented_seq_glnA* > glnA_check.txt
grep -c '>' oriented_seq_acrB* > acrB_check.txt
cd ../

cd BBmap_out_amtB_sbcC
grep -c '>' oriented_seq_amtB* > amtB_check.txt
grep -c '>' oriented_seq_sbcC* > sbcC_check.txt
cd ../

cd BBmap_out_glnA_helY
grep -c '>' oriented_seq_glnA* > glnA_check.txt
grep -c '>' oriented_seq_helY* > helY_check.txt
cd ../

cd BBmap_out_helY_sbcC
grep -c '>' oriented_seq_helY* > helY_check.txt
grep -c '>' oriented_seq_sbcC* > sbcC_check.txt
cd ../



# Run a python script written by Chris McLimans that checks the number of sequences extracted per sample from
# the written grep files above and if no sequence is written then removes the sample from the master
# file list (files.txt) and writes a new list of sample IDs to the _markers.txt file for samples that have
# a sequence.
#
# NOTE: the python script currently requires that the first sample in the check has the format:
# oriented_seq_glnA_M025S2.fa with 3 underscores ( _ ) in the sample ID and before the sample ID.
#

python -u Check_missing.py -i BBmap_out_acrB_glnA/files.txt -p1 BBmap_out_acrB_glnA/glnA_check.txt \
-p2 BBmap_out_acrB_glnA/acrB_check.txt \
-o BBmap_out_acrB_glnA/acrB_glnA_markers.txt


# This if condition checkes if the _markers.txt file is not written then an error occured in the python script
# could be that >1 sequence was extracted for a sample
{
if [ ! -f 'BBmap_out_acrB_glnA/acrB_glnA_markers.txt' ]; then
    echo "acrB_glnA_markers.txt not found!"
    exit 1;
fi
}


python -u Check_missing.py -i BBmap_out_amtB_sbcC/files.txt -p1 BBmap_out_amtB_sbcC/amtB_check.txt \
-p2 BBmap_out_amtB_sbcC/sbcC_check.txt \
-o BBmap_out_amtB_sbcC/amtB_sbcC_markers.txt

{
if [ ! -f 'BBmap_out_amtB_sbcC/amtB_sbcC_markers.txt' ]; then
    echo "amtB_sbcC_markers.txt not found!"
    exit 1;
fi
}



python -u Check_missing.py -i BBmap_out_glnA_helY/files.txt -p1 BBmap_out_glnA_helY/glnA_check.txt \
-p2 BBmap_out_glnA_helY/helY_check.txt \
-o BBmap_out_glnA_helY/glnA_helY_markers.txt

{
if [ ! -f 'BBmap_out_glnA_helY/glnA_helY_markers.txt' ]; then
    echo "glnA_helY_markers.txt not found!"
    exit 1;
fi
}


python -u Check_missing.py -i BBmap_out_helY_sbcC/files.txt -p1 BBmap_out_helY_sbcC/helY_check.txt \
-p2 BBmap_out_helY_sbcC/sbcC_check.txt \
-o BBmap_out_helY_sbcC/helY_sbcC_markers.txt

{
if [ ! -f 'BBmap_out_helY_sbcC/helY_sbcC_markers.txt' ]; then
    echo "helY_sbcC_markers.txt not found!"
    exit 1;
fi
}


# Go into each primer pair directory with the _markers.txt file now written
# Concatenate the individual sequence files together based on the files
# that have a single sequence per the checked _markers.txt file.
# Important: this writes them in order of the _markers.txt file
cd BBmap_out_acrB_glnA
while read line; do cat oriented_seq_acrB_$line >> oriented_acrB_cat.fa; done < acrB_glnA_markers.txt
while read line; do cat oriented_seq_glnA_$line >> oriented_glnA_cat.fa; done < acrB_glnA_markers.txt

cd ../BBmap_out_amtB_sbcC
while read line; do cat oriented_seq_sbcC_$line >> oriented_sbcC_cat.fa; done < amtB_sbcC_markers.txt
while read line; do cat oriented_seq_amtB_$line >> oriented_amtB_cat.fa; done < amtB_sbcC_markers.txt

cd ../BBmap_out_glnA_helY
while read line; do cat oriented_seq_glnA_$line >> oriented_glnA_cat.fa; done < glnA_helY_markers.txt
while read line; do cat oriented_seq_helY_$line >> oriented_helY_cat.fa; done < glnA_helY_markers.txt

cd ../BBmap_out_helY_sbcC
while read line; do cat oriented_seq_sbcC_$line >> oriented_sbcC_cat.fa; done < helY_sbcC_markers.txt
while read line; do cat oriented_seq_helY_$line >> oriented_helY_cat.fa; done < helY_sbcC_markers.txt
cd ../


# This bash script written by Chris McLimans checked the number of sequneces expected
# per the _markers.txt files for a marker pair against the number of sequences
# in the concatenated file for each primer within the pair.
# Need to check the stdout file to confirm all are equal
./check_counts.sh BBmap_out_acrB_glnA/acrB_glnA_markers.txt BBmap_out_acrB_glnA/oriented_acrB_cat.fa
./check_counts.sh BBmap_out_acrB_glnA/acrB_glnA_markers.txt BBmap_out_acrB_glnA/oriented_glnA_cat.fa

./check_counts.sh BBmap_out_amtB_sbcC/amtB_sbcC_markers.txt BBmap_out_amtB_sbcC/oriented_sbcC_cat.fa
./check_counts.sh BBmap_out_amtB_sbcC/amtB_sbcC_markers.txt BBmap_out_amtB_sbcC/oriented_amtB_cat.fa

./check_counts.sh BBmap_out_glnA_helY/glnA_helY_markers.txt BBmap_out_glnA_helY/oriented_glnA_cat.fa
./check_counts.sh BBmap_out_glnA_helY/glnA_helY_markers.txt BBmap_out_glnA_helY/oriented_helY_cat.fa

./check_counts.sh BBmap_out_helY_sbcC/helY_sbcC_markers.txt BBmap_out_helY_sbcC/oriented_sbcC_cat.fa
./check_counts.sh BBmap_out_helY_sbcC/helY_sbcC_markers.txt BBmap_out_helY_sbcC/oriented_helY_cat.fa


### Step 4
# Replace headers with file names
#
# Use a python script written by Chris McLimans to replace the fasta headers in the cat files
# with the sample ID specified in the _marker files for each primer pair.
# The sequence order in the cat files follows the _markers.txt order and thus can be followed to
# replace the headers with sample IDs.
###

./Replace_headers_commandline.py -m BBmap_out_acrB_glnA/acrB_glnA_markers.txt -i BBmap_out_acrB_glnA/oriented_acrB_cat.fa -o BBmap_out_acrB_glnA/header_oriented_acrB_cat.fa
./Replace_headers_commandline.py -m BBmap_out_acrB_glnA/acrB_glnA_markers.txt -i BBmap_out_acrB_glnA/oriented_glnA_cat.fa -o BBmap_out_acrB_glnA/header_oriented_glnA_cat.fa

./Replace_headers_commandline.py -m BBmap_out_amtB_sbcC/amtB_sbcC_markers.txt -i BBmap_out_amtB_sbcC/oriented_sbcC_cat.fa -o BBmap_out_amtB_sbcC/header_oriented_sbcC_cat.fa
./Replace_headers_commandline.py -m BBmap_out_amtB_sbcC/amtB_sbcC_markers.txt -i BBmap_out_amtB_sbcC/oriented_amtB_cat.fa -o BBmap_out_amtB_sbcC/header_oriented_amtB_cat.fa

./Replace_headers_commandline.py -m BBmap_out_glnA_helY/glnA_helY_markers.txt -i BBmap_out_glnA_helY/oriented_glnA_cat.fa -o BBmap_out_glnA_helY/header_oriented_glnA_cat.fa
./Replace_headers_commandline.py -m BBmap_out_glnA_helY/glnA_helY_markers.txt -i BBmap_out_glnA_helY/oriented_helY_cat.fa -o BBmap_out_glnA_helY/header_oriented_helY_cat.fa

./Replace_headers_commandline.py -m BBmap_out_helY_sbcC/helY_sbcC_markers.txt -i BBmap_out_helY_sbcC/oriented_sbcC_cat.fa -o BBmap_out_helY_sbcC/header_oriented_sbcC_cat.fa
./Replace_headers_commandline.py -m BBmap_out_helY_sbcC/helY_sbcC_markers.txt -i BBmap_out_helY_sbcC/oriented_helY_cat.fa -o BBmap_out_helY_sbcC/header_oriented_helY_cat.fa






# Based on OSCER software, unload the BBMap software and then start
# muscle to perform alignments
module unload BBMap
source activate muscle

### Step 5
# Do Alignments
#
# Use muscle to align the cat files with headers replace with sample IDs
# Outputs an alignment for each primer within a primer pair into a new
# align directory within the primer pair directory
#
# Finally, copy the FASconCAT software (or acquire from github) into the
# align directory
###

cd BBmap_out_acrB_glnA/
mkdir align
muscle -in header_oriented_acrB_cat.fa -clwstrict -out align/acrB.aln
muscle -in header_oriented_glnA_cat.fa -clwstrict -out align/glnA.aln
cp ~/script_archive/FASconCAT.pl align/
cd ../

cd BBmap_out_amtB_sbcC/
mkdir align
muscle -in header_oriented_sbcC_cat.fa -clwstrict -out align/sbcC.aln
muscle -in header_oriented_amtB_cat.fa -clwstrict -out align/amtB.aln
cp ~/script_archive/FASconCAT.pl align/
cd ../

cd BBmap_out_glnA_helY/
mkdir align
muscle -in header_oriented_glnA_cat.fa -clwstrict -out align/glnA.aln
muscle -in header_oriented_helY_cat.fa -clwstrict -out align/helY.aln
cp ~/script_archive/FASconCAT.pl align/
cd ../

cd BBmap_out_helY_sbcC/
mkdir align
muscle -in header_oriented_sbcC_cat.fa -clwstrict -out align/sbcC.aln
muscle -in header_oriented_helY_cat.fa -clwstrict -out align/helY.aln
cp ~/script_archive/FASconCAT.pl align/
cd ../


# Close the muscle software
conda deactivate


### Step 6
# FASconCAT super alignments
#
# use FASconCAT to combine both primer pairs into one sequence per sample
###

cd BBmap_out_acrB_glnA/align/
perl FASconCAT.pl -o -o -s
cd ../../

cd BBmap_out_amtB_sbcC/align/
perl FASconCAT.pl -o -o -s
cd ../../

cd BBmap_out_glnA_helY/align/
perl FASconCAT.pl -o -o -s
cd ../../

cd BBmap_out_helY_sbcC/align/
perl FASconCAT.pl -o -o -s
cd ../../


###
# RAxML Trees
###

module load RAxML/8.2.11-intel-2018a-pthreads-avx2

cd BBmap_out_acrB_glnA/align/
raxmlHPC -f a -T 8 -s FcC_supermatrix.fas -o Ma_SC_T_19800800_S464 -n Expanded_acrB_glnA_tree -w scratch/cmclimans/botrys/$1/BBmap_out_acrB_glnA/align/ -m GTRCAT -x 522 -p 22 -# 200
cd ../../

cd BBmap_out_amtB_sbcC/align/
raxmlHPC -f a -T 8 -s FcC_supermatrix.fas -o Ma_SC_T_19800800_S464 -n Expanded_amtB_sbcC_tree -w scratch/cmclimans/botrys/$1/BBmap_out_amtB_sbcC/align/ -m GTRCAT -x 522 -p 22 -# 200
cd ../../

cd BBmap_out_glnA_helY/align/
raxmlHPC -f a -T 8 -s FcC_supermatrix.fas -o Ma_SC_T_19800800_S464 -n Expanded_glnA_helY_tree -w scratch/cmclimans/botrys/$1/BBmap_out_glnA_helY/align/ -m GTRCAT -x 522 -p 22 -# 200
cd ../../

cd BBmap_out_helY_sbcC/align/
raxmlHPC -f a -T 8 -s FcC_supermatrix.fas -o Ma_SC_T_19800800_S464 -n Expanded_helY_sbcC_tree -w scratch/cmclimans/botrys/$1/BBmap_out_helY_sbcC/align/ -m GTRCAT -x 522 -p 22 -# 200
cd ../../
