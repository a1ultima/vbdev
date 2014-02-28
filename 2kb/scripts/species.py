def generate_list(speciesFile='./species_list.txt'):
    print 'loading species...'
    file_in = open(speciesFile)
    speciesList = []
    while True:
        line = file_in.readline()
        if line == "":                    
            break
        print('\t'+line.rstrip())
        line_split = line.split(' ')      # split line to separate comments away from species name
        speciesList.append(line_split[0].rstrip())
    file_in.close()
    return speciesList

if __name__ == "__main__":
    import sys
    generate_list(sys.argv[1])