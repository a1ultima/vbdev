
import os
import subprocess
from datetime import datetime
import speciesManage
import sys

def allSpecies( speciesFile, motifsPath_in, fastaPath_in, dataPath_out, verbosity, threshold ):
    """
    Notes:

    Args:
        speciesFile     = './species_list.txt'
        motifsPath_in   = '../data/meme_data/out/dreme_100bp/sampled_all_hpc/'
        fastaPath_in    = '../data/meme_data/in/random_dreme/'
        dataPath_out    = '../data/meme_data/out/fimo/'
        verbosity       = 2                                                     # 
        threshold       = 0.1                                                   # q-value cut-off for bad matches filter
    """

    species_list    = speciesManage.generate_list(speciesFile) # generates a list of species names corresponding to EnsEMBl MySQL db names

    if not os.path.exists(dataPath_out):
    	os.makedirs(dataPath_out)

    suffix_in_motifs= '_100bp/dreme.html'
    suffix_in_fasta = '_upstream_dremeready_all_simpleMasked_random.fasta'
    suffix_out      = '_100bp'
    
    for species in species_list:
        print(species)
        motifs_in = motifsPath_in+ species + suffix_in_motifs
        fasta_in = fastaPath_in + species + suffix_in_fasta
        fimo_out = dataPath_out + species + suffix_out
        # COMMAND LINE: fimo -oc fimo_100bp/ -thresh 0.1 -verbosity 2 ./dreme_100bp/sampled_all_hpc/anopheles_christyi_100bp/dreme.html ./dreme_100bp/sampled_all_hpc/anopheles_merus_100bp/dreme.html
        fimo = [    'fimo',                           # fimo program
                    '-verbosity',   str(verbosity),     # verbosity, 1:5
                    '-oc',          fimo_out,          # over-write directory with <name> and write in the results
                    '-thresh',      str(threshold),     # significance threshold, q-value
                    motifs_in,                  # positive file: input fasta data
                    fasta_in
                    ]
        all_fimo = [fimo]
        for fimo in all_fimo:
            #start   = subprocess.Popen(['date'], stdout=subprocess.PIPE).communicate()[0]
            subprocess.call(fimo)
            #end     = subprocess.Popen(['date'], stdout=subprocess.PIPE).communicate()[0]
            #tdelta  = datetime.strptime(end.split(' ')[3], '%H:%M:%S') - datetime.strptime(start.split(' ')[3], '%H:%M:%S')
            #print(tdelta)

if __name__ == "__main__":
    try: 
        speciesFile     = sys.argv[1]
    except IndexError:
        speciesFile     = './species_list.txt'
    try: 
        motifsPath_in   = sys.argv[2]
    except IndexError:
        motifsPath_in   = '../data/meme_data/out/dreme_100bp/sampled_all_hpc/'
    try:
        fastaPath_in    = sys.argv[3]
    except IndexError:
        fastaPath_in    = '../data/meme_data/in/random_dreme/'
    try: 
        dataPath_out    = sys.argv[4]
    except IndexError:
        dataPath_out    = '../data/meme_data/out/fimo_100bp/'
    try: 
        verbosity       = int(sys.argv[5])
    except IndexError:
        verbosity       = 2
    try: 
        threshold       = float(sys.argv[6])
    except IndexError:
        threshold       = 0.1
    allSpecies(speciesFile, motifsPath_in, fastaPath_in, dataPath_out, verbosity, threshold)
    print('COMPLETE!')