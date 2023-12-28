#!/bin/bash
#
#SBATCH --partition=debug
#SBATCH --ntasks=1
#SBATCH --mem=1GB
#SBATCH --output=jobname_%J_stdout.txt
#SBATCH --error=jobname_%J_stderr.txt
#SBATCH --time=00:30:00
#SBATCH --job-name=SAMtools
#SBATCH --mail-user=cmclimans@ou.edu
#SBATCH --chdir=/home/cmclimans/
#
#################################################


module load SAMtools
export tmp_PATH='/scratch/cmclimans/prep_analysis/single_geno_genes/reformatted'

#samtools faidx $tmp_PATH/seqs_orthos_single_orthos_B.txt $(cat $tmp_PATH/orthos_single_orthos_B.txt) > $tmp_PATH/reordered_seqs_orthos_single_orthos_B.txt
#samtools faidx $tmp_PATH/seqs_orthos_single_orthos_C.txt $(cat $tmp_PATH/orthos_single_orthos_C.txt) > $tmp_PATH/reordered_seqs_orthos_single_orthos_C.txt
#samtools faidx $tmp_PATH/seqs_orthos_single_orthos_E.txt $(cat $tmp_PATH/orthos_single_orthos_E.txt) > $tmp_PATH/reordered_seqs_orthos_single_orthos_E.txt
#samtools faidx $tmp_PATH/seqs_orthos_single_orthos_F.txt $(cat $tmp_PATH/orthos_single_orthos_F.txt) > $tmp_PATH/reordered_seqs_orthos_single_orthos_F.txt
#samtools faidx $tmp_PATH/seqs_orthos_single_orthos_G.txt $(cat $tmp_PATH/orthos_single_orthos_G.txt) > $tmp_PATH/reordered_seqs_orthos_single_orthos_G.txt
#samtools faidx $tmp_PATH/seqs_orthos_single_orthos_H.txt $(cat $tmp_PATH/orthos_single_orthos_H.txt) > $tmp_PATH/reordered_seqs_orthos_single_orthos_H.txt

samtools faidx $tmp_PATH/seqs_orthos_single_orthos_I.txt $(cat $tmp_PATH/orthos_single_orthos_I.txt) > $tmp_PATH/reordered_seqs_orthos_single_orthos_I.txt
samtools faidx $tmp_PATH/seqs_orthos_single_orthos_J.txt $(cat $tmp_PATH/orthos_single_orthos_J.txt) > $tmp_PATH/reordered_seqs_orthos_single_orthos_J.txt
samtools faidx $tmp_PATH/seqs_orthos_single_orthos_K.txt $(cat $tmp_PATH/orthos_single_orthos_K.txt) > $tmp_PATH/reordered_seqs_orthos_single_orthos_K.txt
samtools faidx $tmp_PATH/seqs_orthos_single_orthos_L.txt $(cat $tmp_PATH/orthos_single_orthos_L.txt) > $tmp_PATH/reordered_seqs_orthos_single_orthos_L.txt
samtools faidx $tmp_PATH/seqs_orthos_single_orthos_M.txt $(cat $tmp_PATH/orthos_single_orthos_M.txt) > $tmp_PATH/reordered_seqs_orthos_single_orthos_M.txt
samtools faidx $tmp_PATH/seqs_orthos_single_orthos_N.txt $(cat $tmp_PATH/orthos_single_orthos_N.txt) > $tmp_PATH/reordered_seqs_orthos_single_orthos_N.txt
samtools faidx $tmp_PATH/seqs_orthos_single_orthos_O.txt $(cat $tmp_PATH/orthos_single_orthos_O.txt) > $tmp_PATH/reordered_seqs_orthos_single_orthos_O.txt
samtools faidx $tmp_PATH/seqs_orthos_single_orthos_P.txt $(cat $tmp_PATH/orthos_single_orthos_P.txt) > $tmp_PATH/reordered_seqs_orthos_single_orthos_P.txt
samtools faidx $tmp_PATH/seqs_orthos_single_orthos_Q.txt $(cat $tmp_PATH/orthos_single_orthos_Q.txt) > $tmp_PATH/reordered_seqs_orthos_single_orthos_Q.txt
samtools faidx $tmp_PATH/seqs_orthos_single_orthos_R.txt $(cat $tmp_PATH/orthos_single_orthos_R.txt) > $tmp_PATH/reordered_seqs_orthos_single_orthos_R.txt
samtools faidx $tmp_PATH/seqs_orthos_single_orthos_cluster2O.txt $(cat $tmp_PATH/orthos_single_orthos_cluster2O.txt) > $tmp_PATH/reordered_seqs_orthos_single_orthos_cluster2O.txt
samtools faidx $tmp_PATH/seqs_orthos_single_orthos_cluster18.txt $(cat $tmp_PATH/orthos_single_orthos_cluster18.txt) > $tmp_PATH/reordered_seqs_orthos_single_orthos_cluster18.txt
samtools faidx $tmp_PATH/seqs_orthos_single_orthos_cluster19.txt $(cat $tmp_PATH/orthos_single_orthos_cluster19.txt) > $tmp_PATH/reordered_seqs_orthos_single_orthos_cluster19.txt
samtools faidx $tmp_PATH/seqs_orthos_single_orthos_cluster21.txt $(cat $tmp_PATH/orthos_single_orthos_cluster21.txt) > $tmp_PATH/reordered_seqs_orthos_single_orthos_cluster21.txt
samtools faidx $tmp_PATH/seqs_orthos_single_orthos_clusterO5.txt $(cat $tmp_PATH/orthos_single_orthos_clusterO5.txt) > $tmp_PATH/reordered_seqs_orthos_single_orthos_clusterO5.txt
samtools faidx $tmp_PATH/seqs_orthos_single_orthos_clusterO7.txt $(cat $tmp_PATH/orthos_single_orthos_clusterO7.txt) > $tmp_PATH/reordered_seqs_orthos_single_orthos_clusterO7.txt
samtools faidx $tmp_PATH/seqs_orthos_single_orthos_clusterO8.txt $(cat $tmp_PATH/orthos_single_orthos_clusterO8.txt) > $tmp_PATH/reordered_seqs_orthos_single_orthos_clusterO8.txt
