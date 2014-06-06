#!/bin/bash

# Before running:
#   - scp the hpc data into         : /data/meme_data/out/dreme_100bp/hpc_scp then check it
#   - copy over to sampled_all_hpc  : cp /data/meme_data/out/dreme_100bp/hpc_scp/* /data/meme_data/out/dreme_100bp/sampled_all_hpc/evalue_0.05/
#   - run motifPmw_toFasta.py       : >>>import motifPmw_toFasta >>> cool = motifPmw_toFasta.allSpecies( speciesFile='./species_list.txt', dataPath_in='../data/meme_data/out/dreme_100bp/sampled_all_hpc/evalue_0.05/', method='TRANSFAC' )
#   - copy all_100bp/dreme.fasta    : from stamp_data/ : cp ../meme_data/out/dreme_100bp/sampled_all_hpc/evalue_0.05/all_100bp/dreme.fasta ./in/dreme_100bp_e0.05/
#   - run stamp
#

# 0. GOTO   : this can be altered
cd ~/0VB/2kb/scripts/ &&

# 1. COPY   : dreme output HPC --> VB-DEV
scp -r ab108@login.cx1.hpc.ic.ac.uk:/home/ab108/work/dreme_testing_e0.05/out/*_100bp ../data/meme_data/out/dreme_100bp/hpc_scp/ &&
cp ../data/meme_data/out/dreme_100bp/hpc_scp/*_100bp ../data/meme_data/out/dreme_100bp/sampled_all_hpc/evalue_0.05/ &&

# 2. UNITE  : unite all species TRANSFAC
python motifPmw_toFasta.py ./species_list.txt ../data/meme_data/out/dreme_100bp/sampled_all_hpc/evalue_0.05/ TRANSFAC &&

# 3. STAMP  : prep & run


cp ../data/meme_data/out/dreme_100bp/sampled_all_hpc/evalue_0.05/all_100bp/dreme.fasta ../data/stamp_data/in/dreme_100bp_e0.05/ &&

# DIRECT    : but error prone, maybe stamp can't interpret the . in .05 ??
MATCH=../data/stamp_data/in/common/jaspar.motifs
TF=../data/stamp_data/in/dreme_100bp_e0.05/dreme.fasta
SD=../data/stamp_data/in/common/ScoreDists/JaspRand_SSD_SWU.scores
OUT=../data/stamp_data/out/dreme_100bp_e0.05/SWU_SSD/out

stamp -tf $TF -sd $SD -match $MATCH -cc SSD -align SWU -out $OUT

# INDIRECT



#echo COMPLETE

