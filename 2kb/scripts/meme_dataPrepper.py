import subprocess

def memedata_prep( filename_in, filename_out ): 

    # Need to ensure data is ready for MEME
    file_in  = open(filename_in)        # Read file in order to establish filters that will later format the data ready for MEME
    file_out = open(filename_out,'w')
    headers = {}

    print '\tEstablishing filters...'
    while True:
        header  = file_in.readline()
        sequence= file_in.readline()
        if header == "":                    # break when finished
            break
        header_split    = header.split('\t')
        geneId          = header_split[0].replace('>','')
        utrCoord        = header_split[1].replace('UtrCoord:','').split(':')
        utrLength       = header_split[2].replace('UtrLength:','')
        utrType         = header_split[3].replace('UtrType:','').rstrip()
        chromosome      = utrCoord[0]

        startEnd        = utrCoord[1].rpartition('-')   # rpartition ensures: (-800-20) -> '-800','-','20' 
        start           = startEnd[0]
        end             = startEnd[2]

        strand          = utrCoord[2]
        headers[geneId] = {'utrCoord':{'chromosome':chromosome,'start':start,'end':end,'strand':strand},'utrLength':utrLength,'utrType':utrType}
    file_in.close()

    # Pulley samples:       shared samples with genes in the same direction 
    #
    #   ---S--- + ---G-------->
    #   ---S--- + ---G--->
    #   ---S--- + ---G------------>     ...only one of these samples are kept, rest are thrown...
    # 
    #   dealt with by making a dictionary that ensures a 1 to 1 key-value relationshipof sample-to-geneId
    #
    #   during the next .readline() loop utrcoord-to-geneid value will be checked against the current geneId in current iteration
    #
    geneIds = headers.keys()
    one_coord_per_geneId = {}
    for geneId in geneIds:


        utrCoord = headers[geneId]['utrCoord']['chromosome']+':'+headers[geneId]['utrCoord']['start']+'-'+headers[geneId]['utrCoord']['end']+':'+headers[geneId]['utrCoord']['strand']


        one_coord_per_geneId[utrCoord] = geneId   # To ensure only one geneId = one utrCoord, I purposefully allow values to be overidden for a given key 

    # Tug'o'war samples:    shared samples with genes in the opposite direction 
    #
    #   s1  <-------G--- + ---S---
    #   s2                 ---S--- + ---S------->
    #
    #   dealt with by making a set of unique utrcoords stripped of their strand property, now upon encountering a match forward samples and reverse samples take their half of the tug-o-war sequence
    #
    coords          = [ ':'.join(i.split(':')[0:-1]) for i in one_coord_per_geneId.keys() ]     # remove strand property
    coords_shared   = set([i for i in coords if coords.count(i)>1])

    print '\tGenerating meme-ready data...'
    file_in = open(filename_in) # Open file for formatted parsing
    collection = []
    coords_shared_geneIds = []

    print ('\t\tFiltering shared samples...')
    while True:    
        header  = file_in.readline()
        sequence= file_in.readline().rstrip()
        if header == "": # break when finished
            break
    # MEME filter:  =>  Is this sample seq too short for meme? 
        # YES:
        if len(sequence)<10:
            if sequence=='\n':
                collection.append(header)
            continue
        # NO:
        else:
            geneId  = header.split('\t')[0].replace('>','')
            utrCoord= header.split('\t')[1].replace('UtrCoord:','').split(':')
        # PULLEY filter:            =>  Is this GeneId unique? or permitted as a winning pulley sample?

            #'supercont1.56:-1973-27:1'

            if one_coord_per_geneId[':'.join(utrCoord)] == geneId:
                utrCoord_noStrand   = ':'.join(utrCoord[0:1+1])
            # TUG'O'WAR filter:     =>  Is this sample in a tug'o'war? 
                if utrCoord_noStrand in coords_shared:
                    #print(geneId)
                    #utrCoord_onlyStrand = utrCoord[2]
                    coords_shared_geneIds.append(header.split('\t')[0].replace('>',''))
                    sequence = sequence[(len(sequence)/2):] # this works, since the rightermost end of the seq = gene, so we must sample from the end inwards
                # OKAY:             => then upon passing the filters, we can write the sample :)
                file_out.write(header+sequence+'\n')
    file_in.close()
    file_out.close()    

def memedata_dust ( filename_in, cutoff=20 ): 

    # DUSTMASKER: low complexity region masker
    print('\t\tMasking low complexity regions...')
    file_dust_name  = filename_in.split('.')[0]+'_dusted'+'.'+filename_in.split('.')[1]
    file_dust       = open(file_dust_name,'w')
    subprocess.call(['dust',filename_in,'10'],stdout=file_dust) # TODO: what value do I use for 'cut-off', here given as 10
    file_dust.close()


def allSpecies_memedata ( species_list, dataPath_in='home/ab108/0VB/2kb/data/sample_seqs/fasta/masking/', dataPath_out = 'home/ab108/0VB/2kb/data/meme_data/in/'): 
    """
    Notes:
        prepares the 2kb upstream data for MEME to start processing for every species
    Args:
        species_list        = ['Anopheles gambiae', 'Aedes aegypti']    #list: species names to sample from
    """

    import os

    for species in species_list:
        print(species)

        species = species.lower().replace(' ','_')
        home = os.path.expanduser('~')+'/'

        os.chdir(home)

        dataPath_in = os.path.expanduser(dataPath_in.replace(home[1:],''))
        dataPath_out= os.path.expanduser(dataPath_out.replace(home[1:],''))

        filename_in = dataPath_in + species + '_upstream.fasta'
        filename_out= dataPath_out+ species + '_upstream_memeready_all.fasta'

        memedata_prep(filename_in,filename_out)
        memedata_dust(filename_out)


#-------------------------------------------------------------------------------------------------------------------------------------
# RUN
#-------------------------------------------------------------------------------------------------------------------------------------
species_list = [    'Aedes aegypti',
                'Anopheles albimanus',
                'Anopheles arabiensis',
                'Anopheles atroparvus',
                'Anopheles christyi',
                'Anopheles culicifacies',
                'Anopheles darlingi',
                'Anopheles dirus',
                'Anopheles epiroticus',
                'Anopheles farauti',
                'Anopheles funestus',
                'Anopheles gambiae',
                'Anopheles maculatus',
                'Anopheles melas',
                'Anopheles merus',
                'Anopheles minimus',
                'Anopheles quadriannulatus',
                'Anopheles sinensis'
                #,'Anopheles stephensiI',           # this is broken..?
                #'Anopheles stephensi'              # ''
                ]

""" NOTICE:


"""

allSpecies_memedata (species_list, dataPath_in='home/ab108/0VB/2kb/data/sample_seqs/fasta/masking/', dataPath_out = 'home/ab108/0VB/2kb/data/meme_data/in/')

# # TEST MODE:
# # filename_in    = '../data/meme_data/test/anopheles_gambiae_upstream_memeready_100.fasta'
# # filename_out   = '../data/meme_data/test/anopheles_gambiae_upstream_memeready_independent.fasta'

# # REAL MODE:
# filename_in    = '../data/meme_data/test/anopheles_gambiae_upstream.fasta'
# filename_out   = '../data/meme_data/test/anopheles_gambiae_upstream_memeready_all.fasta'

# memedata_prep(filename_in,filename_out)
# memedata_dust(filename_out,cutoff=10)

print('COMPLETE!')