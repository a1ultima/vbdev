import subprocess
#from datetime import datetime
import speciesManage
import sys

def allSpecies( speciesFile, dataPath_in, dataPath_out, resultFormat, verbosity, eValue, nMotifs ):
    """
    Notes:

    Args:
        speciesFile     = './species_list.txt'
        dataPath_in     = '../data/meme_data/in/random_dreme/'
        dataPath_out    = '../data/meme_data/out/'
        verbosity       = 2
        eValue          = 0.0001
        nMotifs         = 100
        resultFormat    = 'png'
        hpc             = True                                  # is this run via HPC?
    """

    suffix_in       = '_upstream_dremeready_all_simpleMasked_random.fasta'
    suffix_out      = '_100bp'
    species_list    = speciesManage.generate_list(speciesFile) # generates a list of species names corresponding to EnsEMBl MySQL db name

    for species in species_list:

        filename_in = dataPath_in + species + suffix_in

        dreme_out   = dataPath_out+species+suffix_out
        # dreme -png -v 1 -oc . -t 18000 -p anopheles_gambiae_upstream_memeready_all_simpleMasked.fasta -e 0.05 -m 100 -dfile description
        dreme = [   'dreme',                    # DREME program
                    '-'+resultFormat,           # result format
                    '-v', str(verbosity),       # verbosity, 1:5
                    '-oc', dreme_out,           # over-write directory with <name> and write in the results
                    #'-t','18000',              # time-limit
                    '-p', filename_in,          # positive file: input fasta data
                    #'-n', ... ,                # negative file: the null model of fasta sequences, i.e. dinucleotide shuffling 
                    '-e', str(eValue),    # E-value cut-off, ignores the list of motifs with < -e from running the long analysis on
                    '-m', str(nMotifs),         # n_motifs to search for, runtime scales linearly w/ acceleration
                    ]

        # if hpc==True:
        #     dreme = ['qsub']+dreme

        all_dreme = [dreme]

        for dreme in all_dreme:
            #start   = subprocess.Popen(['date'], stdout=subprocess.PIPE).communicate()[0]
            subprocess.call(dreme)
            #end     = subprocess.Popen(['date'], stdout=subprocess.PIPE).communicate()[0]
            #tdelta  = datetime.strptime(end.split(' ')[3], '%H:%M:%S') - datetime.strptime(start.split(' ')[3], '%H:%M:%S')
            #print(tdelta)

if __name__ == "__main__":
    try: 
        speciesFile = sys.argv[1]
    except IndexError:
        speciesFile = './species_list.txt'
    try: 
        dataPath_in = sys.argv[2]
    except IndexError:
        dataPath_in = '../data/meme_data/in/random_dreme/'
    try: 
        dataPath_out = sys.argv[3]
    except IndexError:
        dataPath_out = '../data/meme_data/out/'
    try: 
        resultFormat = sys.argv[4]
    except IndexError:
        resultFormat = 'png'
    try: 
        verbosity = sys.argv[5]
    except IndexError:
        verbosity = 2
    try: 
        eValue = sys.argv[6]
    except IndexError:
        eValue = 0.0001
    try: 
        nMotifs = sys.argv[7]
    except IndexError:
        nMotifs = 100
    # try: 
    #     hpc = bool(sys.argv[8])
    # except IndexError:
    #     hpc = True
    allSpecies(speciesFile, dataPath_in, dataPath_out, resultFormat, verbosity, eValue, nMotifs)
    print('COMPLETE!')