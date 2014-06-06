#!/bin/bash

# Before running:
#   - scp the hpc data into         : /data/meme_data/out/dreme_100bp/hpc_scp then check it
#   - copy over to sampled_all_hpc  : cp /data/meme_data/out/dreme_100bp/hpc_scp/* /data/meme_data/out/dreme_100bp/sampled_all_hpc/evalue_0.05/
#   - run motifPmw_toFasta.py       : >>>import motifPmw_toFasta >>> cool = motifPmw_toFasta.allSpecies( speciesFile='./species_list.txt', dataPath_in='../data/meme_data/out/dreme_100bp/sampled_all_hpc/evalue_0.05/', method='TRANSFAC' )
#   - copy all_100bp/dreme.fasta    : from stamp_data/ : cp ../meme_data/out/dreme_100bp/sampled_all_hpc/evalue_0.05/all_100bp/dreme.fasta ./in/dreme_100bp_e0.05/
#   - run stamp
#

# Backup previous files
shopt -s extglob    # enable extended regex

DATE=05-25
ECUT=0.05
KBUP=2000

#----------------------------------------------------------------------

# Replace input data for dreme_hpc: ../../hpc/dreme/dreme_testing_e0.05/data/*
#
rm -r ../hpc/dreme/dreme_testing_e$ECUT\_b$KBUP/data/*
# 
# Unpack a backup
#
cp -r ../data/meme_data/in/backup-05-12_1500bp_promoters/* ../data/meme_data/in/
#
# Copy dreme input data to the HPC-to-go then shp it off
#
cp -r ../data/meme_data/in/random_dreme/* ../hpc/dreme/dreme_testing_e$ECUT\_b$KBUP/data/
scp -r ../hpc/dreme/dreme_testing_e$ECUT\_b$KBUP/ ab108@login.cx1.hpc.ic.ac.uk:/home/ab108/work/

# Backup Core pipeline data: /meme_data/in
#
# mkdir ../data/meme_data/in/backup-$DATE/
# mv ../data/meme_data/in/!(backup*) ../data/meme_data/in/backup-$DATE/

# Restore Core pipeline data: /meme_data/in
#
# mv ../data/meme_data/in/backup-$DATE/* ../data/meme_data/in/
#
# OR
#
# mv ../data/meme_data/in/backup-$DATE/* ../data/meme_data/in/

#----------------------------------------------------------------------
# Execute the following pipeline when dreme_hpc is finished
#----------------------------------------------------------------------

# BACKUPS
mkdir ../data/stamp_data/out/dreme_100bp_e$ECUT/SWU_SSD/backup-$DATE # create backupfile                                            # ../data/stamp_data/out/dreme_100bp_e$ECUT/SWU_SSD

mv ../data/stamp_data/out/dreme_100bp_e$ECUT/SWU_SSD/!(backup*) ../data/stamp_data/out/dreme_100bp_e$ECUT/SWU_SSD/backup-$DATE/ #move all files except for backup into new backup

# Recover backup ^
# cp -r ../data/stamp_data/out/dreme_100bp_e$ECUT/SWU_SSD/backup-$DATE/* ../data/stamp_data/out/dreme_100bp_e$ECUT/SWU_SSD/ #move all files except for backup into new backup

mkdir ../data/meme_data/out/dreme_100bp/hpc_scp/backup-$DATE # ../data/meme_data/out/dreme_100bp/hpc_scp

mv ../data/meme_data/out/dreme_100bp/hpc_scp/!(backup*) ../data/meme_data/out/dreme_100bp/hpc_scp/backup-$DATE/

mkdir ../data/meme_data/out/dreme_100bp/sampled_all_hpc/evalue_$ECUT/backup-$DATE                                                   # ../data/meme_data/out/dreme_100bp/sampled_all_hpc/evalue_$ECUT/backup-$DATE
mv ../data/meme_data/out/dreme_100bp/sampled_all_hpc/evalue_$ECUT/!(backup*) ../data/meme_data/out/dreme_100bp/sampled_all_hpc/evalue_$ECUT/backup-$DATE/

mkdir ../data/stamp_data/in/dreme_100bp_e$ECUT/backup-$DATE                                                                         # ../data/stamp_data/in/dreme_100bp_e$ECUT/backup-$DATE
mv ../data/stamp_data/in/dreme_100bp_e$ECUT/!(backup*) ../data/stamp_data/in/dreme_100bp_e$ECUT/backup-$DATE/

# 0. GOTO   : this can be altered
cd ~/0VB/2kb/scripts/ 

# 1. COPY FROM HPC  : dreme output HPC --> VB-DEV
scp -r ab108@login.cx1.hpc.ic.ac.uk:/home/ab108/work/dreme_testing_e$ECUT\_b$KBUP/out/*_100bp ../data/meme_data/out/dreme_100bp/hpc_scp/

# everything
#scp -r ab108@login.cx1.hpc.ic.ac.uk:/home/ab108/work/dreme_testing_e$ECUT/ .
#scp -r ab108@login.cx1.hpc.ic.ac.uk:/home/ab108/work/dreme_testing_e$ECUT/data/anopheles_gambiae_upstream_dremeready_all_simpleMasked_random.fasta ../data/meme_data/out/dreme_100bp/hpc_scp/backup/

cp -r ../data/meme_data/out/dreme_100bp/hpc_scp/*_100bp ../data/meme_data/out/dreme_100bp/sampled_all_hpc/evalue_$ECUT/

# 2. UNITE  : unite all species TRANSFAC
python motifPmw_toFasta.py ./species_list.txt ../data/meme_data/out/dreme_100bp/sampled_all_hpc/evalue_$ECUT/ TRANSFAC

# 3. STAMP  : prep & run
cp -r ../data/meme_data/out/dreme_100bp/sampled_all_hpc/evalue_$ECUT/all_100bp/dreme.fasta ../data/stamp_data/in/dreme_100bp_e$ECUT/
MATCH=../data/stamp_data/in/common/jaspar.motifs 
TF=../data/stamp_data/in/dreme_100bp_e$ECUT/dreme.fasta 
SD=../data/stamp_data/in/common/ScoreDists/JaspRand_SSD_SWU.scores 
OUT=../data/stamp_data/out/dreme_100bp_e$ECUT/SWU_SSD/out 
stamp -tf $TF -sd $SD -match $MATCH -cc SSD -align SWU -out $OUT 

# INDIRECT
echo COMPLETE

#--------------------------------------------------------------------
#  BACKUP RECOVERY      <=  use to restore backups for clsuter-analysis pipelines
#--------------------------------------------------------------------

# Recover backup:   dreme_out 
cp -r ../data/stamp_data/in/dreme_100bp_e$ECUT/backup-$DATE/* ../data/stamp_data/in/dreme_100bp_e$ECUT/ # moves contents from a backup to one level up given by parameters: ECUT, DATE

# Recover backup:   stamp_out:  this directory has much of the cluster analysis outputs
cp -r ../data/stamp_data/out/dreme_100bp_e$ECUT/SWU_SSD/backup-$DATE/* ../data/stamp_data/out/dreme_100bp_e$ECUT/SWU_SSD/ #move all files except for backup into new backup



# Recover TREE data for --> cutCollapseStamp 
cp -r ../data/stamp_data/out/dreme_100bp_e0.5/SWU_SSD/backup-05-13/* ../data/stamp_data/out/dreme_100bp_e0.5/SWU_SSD/

