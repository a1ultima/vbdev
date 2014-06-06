import subprocess
from datetime import datetime
import speciesManage
import sys
import itertools

def allSpecies( speciesFile, dataPath_in, dataPath_out, verbosity, threshold ):
    """
    Notes:

    Args:
        speciesFile     = './species_list.txt'
        dataPath_in     = '../data/meme_data/out/dreme_100bp/sampled_all_hpc/'
        dataPath_out    = '../data/meme_data/out/tomtom_100bp/'
        verbosity       = 2                                                     # 
        threshold       = 0.1                                                   # q-value cut-off for bad matches filter
    """
    suffix_in       = '_100bp/dreme.html'
    suffix_out      = '_100bp'
    species_list    = speciesManage.generate_list(speciesFile) # generates a list of species names corresponding to EnsEMBl MySQL db names
    speciesPairs    = list(itertools.product(species_list,species_list))
    for speciesPair in speciesPairs:
        print(speciesPair)
        species_query       = speciesPair[0]
        species_database    = speciesPair[1]
        filename_in_query   = dataPath_in + species_query    + suffix_in
        filename_in_database= dataPath_in + species_database + suffix_in
        tomtom_out          = dataPath_out+species_query+'_vs_'+species_database+suffix_out
        # COMMAND LINE: tomtom -oc tomtom_100bp/ -thresh 0.1 -verbosity 2 ./dreme_100bp/sampled_all_hpc/anopheles_christyi_100bp/dreme.html ./dreme_100bp/sampled_all_hpc/anopheles_merus_100bp/dreme.html
        tomtom = [  'tomtom',                           # Tomtom program
                    '-verbosity',   str(verbosity),     # verbosity, 1:5
                    '-oc',          tomtom_out,          # over-write directory with <name> and write in the results
                    '-thresh',      str(threshold),     # significance threshold, q-value
                    filename_in_query,                  # positive file: input fasta data
                    filename_in_database 
                    ]
        all_tomtom = [tomtom]
        for tomtom in all_tomtom:
            start   = subprocess.Popen(['date'], stdout=subprocess.PIPE).communicate()[0]
            subprocess.call(tomtom)
            end     = subprocess.Popen(['date'], stdout=subprocess.PIPE).communicate()[0]
            tdelta  = datetime.strptime(end.split(' ')[3], '%H:%M:%S') - datetime.strptime(start.split(' ')[3], '%H:%M:%S')
            print(tdelta)

if __name__ == "__main__":
    try: 
        speciesFile = sys.argv[1]
    except IndexError:
        speciesFile = './species_list.txt'
    try: 
        dataPath_in = sys.argv[2]
    except IndexError:
        dataPath_in = '../data/meme_data/out/dreme_100bp/sampled_all_hpc/'
    try: 
        dataPath_out = sys.argv[3]
    except IndexError:
        dataPath_out = '../data/meme_data/out/tomtom_100bp/'
    try: 
        verbosity = sys.argv[4]
    except IndexError:
        verbosity = 2
    try: 
        threshold = sys.argv[5]
    except IndexError:
        threshold = 0.1
    allSpecies(speciesFile, dataPath_in, dataPath_out, verbosity, threshold)
    print('COMPLETE!')