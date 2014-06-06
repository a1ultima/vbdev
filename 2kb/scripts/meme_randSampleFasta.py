import pickle
import os.path
import random
import sys

def allSpecies( speciesFile = './species_list.txt', dataPath_in = '../data/meme_data/in/', dataPath_out='../data/meme_data/in/randomFasta/', n_replicates=5, n_seqs=2000 ):

    import speciesManage # species.generate_list() method creates a python list with names of species corresponding to MySQL EnsEMBL databases
    species_list = speciesManage.generate_list(speciesFile) # generates a list of species names corresponding to EnsEMBl MySQL db names

    for species in species_list:
        print(species)
        suffix_in   = '_upstream_memeready_all_simpleMasked.fasta'
        filename_in = dataPath_in+species+suffix_in
        file_in     = open(filename_in)
        headers     = {}
        seqs        = {}
        
        dataPath_species = dataPath_out+species+'/'

        if not os.path.exists(dataPath_species):
            os.makedirs(dataPath_species) # New dir for species
        else:
            print('\tSamples already AVAILABLE for: '+species)
            continue

        # IMPORT POPULATION FASTA
        print '\tImporting FASTA...'
        while True:
            header  = file_in.readline()
            sequence= file_in.readline()
            if header == '':                    # break when finished
                break
            header_split    = header.split('\t')
            geneId          = header_split[0].replace('>','')
            headers[geneId] = header
            seqs[geneId]    = sequence
        file_in.close()

        # RANDOM SAMPLING LIST
        if os.path.isfile(dataPath_species+'randomList_geneIds.p'):
            print('\tRandom Fasta: AVAILABLE -> Importing...')
            randList_geneIds= pickle.load(open(dataPath_species+'randomList_geneIds.p','rb'))
        else:
            print('\tRandom Fasta: UNAVAILABLE -> Generating...')
            randList        = [[random.randrange(0,len(seqs.keys())-1) for seq in range(0,n_seqs)] for replicate in range(0,n_replicates)]
            randList_geneIds= [[headers.keys()[seq] for seq in randList[replicate]] for replicate in range(0,len(randList))]
            pickle.dump(randList_geneIds,open(dataPath_species+'randomList_geneIds.p','wb'))

        # GENERATE SAMPLE FASTA
        for i,replicate in enumerate(randList_geneIds):
            suffix_out  = suffix_in.replace('_all_','_'+str(n_seqs)+'_'+str(i)+'_')
            #                                ^repl -->    ^#seqs   ^replicate id
            filename_out= dataPath_species+species+suffix_out
            file_out    = open(filename_out,'w')
            for seq in replicate:
                h = headers[seq]
                s = seqs[seq]
                file_out.write(h+s)
            file_out.close()

#-------------------------------------------------------------------------------------------------------------------------------------
# RUN
#-------------------------------------------------------------------------------------------------------------------------------------

# """ NOTICE

# """

# if __name__ == "__main__":
#     try: 
#         speciesFile = sys.argv[1]
#     except IndexError:
#         speciesFile = './species_list.txt',
#     try:
#         dataPath_in = sys.argv[2]
#     except IndexError:
#         dataPath_in = '../data/meme_data/in/'
#     try:
#         dataPath_out = sys.argv[3]
#     except IndexError:
#         dataPath_out = '../data/meme_data/in/randomFasta/'
#     try:
#         n_replicates = int(sys.argv[4])
#     except IndexError:
#         n_replicates = 5
#     try:
#         n_seqs = sys.argv[5]
#     except IndexError:
#         n_seqs = 2000         
#     allSpecies( speciesFile, dataPath_in, dataPath_out, n_replicates, n_seqs )
#     print('COMPLETE!')

