def generate_list(speciesFile='./species_list.txt'):
    """
    Notes:
        Generate a python list of species names which correspond to those in the EnsEMBL MySQL database

    Args:

    Dependencies:
        -->    
            updownstream.py 
            mene_dataPrepper.py 
            meme_bgfileGen.py 
        <--     
            species_list.txt        # the config file for species names                                 
    """
    print 'loading species...'
    file_in = open(speciesFile)
    speciesList = []
    while True:
        line = file_in.readline()
        if line == "":                    
            break
        if line[0] == '#':
            continue
        print('\t'+line.rstrip())
        line_split = line.split(' ')      # split line to separate comments away from species name
        speciesList.append(line_split[0].rstrip())
    file_in.close()
    return speciesList

if __name__ == "__main__":
    import sys
    try:
        generate_list(sys.argv[1])
    except IndexError:
        generate_list('./species_list.txt')