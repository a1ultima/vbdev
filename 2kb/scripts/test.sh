
# DIRECT    : but error prone, maybe stamp can't interpret the . in .05 ??
MATCH=../data/stamp_data/in/common/jaspar.motifs
TF=../data/stamp_data/in/dreme_100bp_e0.05/dreme.fasta
SD=../data/stamp_data/in/common/ScoreDists/JaspRand_SSD_SWU.scores
OUT=../data/stamp_data/out/dreme_100bp_e0.05/SWU_SSD/out

stamp -tf $TF -sd $SD -match $MATCH -cc SSD -align SWU -out $OUT