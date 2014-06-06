#!/bin/bash
#PBS -l select=2:ncpus=16:mem=30000mb
#PBS -l walltime=72:00:00
#PBS -N MEME32

#
# to run it simply type
#
# qsub meme-master.sh
#
# qsub -v SPECIES=anopehels_gambiae meme-master.sh
#

module add intel-suite
module add mpi
module add meme/4.9.1


SPECIES_LIST = 

drosophila_melanogaster
aedes_aegypti
anopheles_albimanus
anopheles_arabiensis
anopheles_atroparvus
anopheles_christyi
anopheles_culicifacies
anopheles_darlingi
anopheles_dirus
anopheles_epiroticus
anopheles_farauti
anopheles_funestus
anopheles_gambiae
anopheles_maculatus
anopheles_melas
anopheles_merus
anopheles_minimus
anopheles_quadriannulatus
anopheles_sinensis
#anopheles_stephI
anopheles_stephensi

SUFFIX_IN=_upstream_dremeready_all_simpleMasked_random.fasta
SUFFIX_OUT=_upstream_dremeready_all_simpleMasked_random.meme.txt

for SPECIES in SPECIES_LIST:

    #SPECIES=anopheles_gambiae
    
    # $PBS_O_WORKDIR is the directory you were in when you did "qsub meme-master.sh"

    # copy all big files from there to this directory
    cp $PBS_O_WORKDIR/data/$SPECIES$SUFFIX_IN infile.fasta

    # copy bgfile too
    #cp $PBS_O_WORKDIR/data/$SPECIES.bg2 bfile.bg2

    # if you run a command with pbsexec -grace 15
    # then it will stop 15 mins before end of job
    # so you could copy checkpoint files back to NFS

    # see if mpiexec is faster
    ls -l

    #mpiexec meme_p infile.fasta -dna -revcomp -mod anr -minsites 200 -maxsites 1000 -nmotifs 100 -w 12 -text -maxsize 50000000 -evt 0.001 -p 32 > outfile.txt
    dreme -png -v 2 -oc dreme_out -p infile.fasta -e 0.0001 -m 1000 

    # copy the results back
    cp outfile.txt $PBS_O_WORKDIR/out/$SPECIES$SUFFIX_OUT

