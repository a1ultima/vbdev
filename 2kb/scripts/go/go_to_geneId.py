# IMPORTS 
import sys
if not '/home/ab108/0VB/2kb/scripts/' in sys.path: # add the 
    sys.path.append('/home/ab108/0VB/2kb/scripts/')
import speciesManage # species.generate_list() method creates a python list with names of species corresponding to MySQL EnsEMBL databases

def allSpecies( speciesFile, dataPath_in, dataPath_out, suffix_in='.txt', suffix_out='_go2geneid.txt'):
    """
    Notes:
        Takes geneId:GoTerm data from Biomart (VectorBase) and generates GoTerm:geneId data

    Args:
        speciesFile = '../species_list.txt'
        dataPath_in = '../../data/go/'
        species     = 'anopheles_gambiae'
    """

    species_list = speciesManage.generate_list(speciesFile) # generates a list of species names corresponding to EnsEMBl MySQL db names

    for species in species_list:

        filename_in = dataPath_in+species+suffix_in
        filename_out= dataPath_in+species+suffix_out
        file_in     = open(filename_in)
        file_out    = open(filename_out,'w')

        # GO:GENEID DICT
        go_to_geneId= {}
        while True:
            line = file_in.readline()
            if line == "":                    # break when finished
                break
            line_split  = line.split('\t')    # split line where '\t' occurs
            geneId      = line_split[0].rstrip()
            goTerm      = line_split[1].rstrip()
            if not goTerm == '':
                if go_to_geneId.has_key(goTerm):
                    go_to_geneId[goTerm].append(geneId)
                else:
                    go_to_geneId[goTerm]=[geneId]
        file_in.close()

        # PARSE FILE
        for goTerm in go_to_geneId.keys():
            file_out.write(goTerm+' '+' '.join(go_to_geneId[goTerm])+'\n')
        file_out.close()

if __name__ == "__main__":
    try: 
        speciesFile = sys.argv[1]
    except IndexError:
        speciesFile = '../species_list.txt'
    try: 
        dataPath_in = sys.argv[2]
    except IndexError:
        dataPath_in = '../data/go/'
    try: 
        dataPath_out= sys.argv[3]
    except IndexError:
        dataPath_out= '../data/go/'
    allSpecies(speciesFile,dataPath_in,dataPath_out)
    print('COMPLETE!')