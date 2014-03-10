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

SPECIES=anopheles_gambiae
SUFFIX=_upstream_memeready_all_dusted

# $PBS_O_WORKDIR is the directory you were in when you did "qsub meme-master.sh"

# copy all big files from there to this directory
cp $PBS_O_WORKDIR/data/$SPECIES$SUFFIX.fasta infile.fasta

# copy bgfile too
cp $PBS_O_WORKDIR/data/$SPECIES.bg2 bfile.bg2

# if you run a command with pbsexec -grace 15
# then it will stop 15 mins before end of job
# so you could copy checkpoint files back to NFS

# see if mpiexec is faster
ls -l

mpiexec meme_p infile.fasta -dna -revcomp -mod anr -minsites 200 -maxsites 1000 -nmotifs 100 -w 12 -text -maxsize 50000000 -evt 0.001 -p 32 > outfile.txt
#meme infile.fasta -p 16 -dna -mod anr -nmotifs 10 -bfile bfile.bg2 -revcomp -text -maxsize 50000000 -evt 0.0 -oc out/

# copy the results back
cp outfile.txt $PBS_O_WORKDIR/out/$SPECIES$SUFFIX.meme.txt

