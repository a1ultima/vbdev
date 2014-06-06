#!/bin/bash

for SPECIES in drosophila_melanogaster aedes_aegypti anopheles_albimanus anopheles_arabiensis anopheles_atroparvus anopheles_christyi anopheles_culicifacies anopheles_darlingi anopheles_dirus anopheles_epiroticus anopheles_farauti anopheles_funestus anopheles_gambiae anopheles_maculatus anopheles_melas anopheles_merus anopheles_minimus anopheles_quadriannulatus anopheles_sinensis anopheles_stephensi

do 
	echo $SPECIES
    # qsub -v SPECIES=$SPECIES ./dreme_hpc.sh
done