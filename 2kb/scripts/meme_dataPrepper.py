import subprocess
import re 
import os

def memedata_prep( filename_in, filename_out ): 
    """
    Notes:
        Preparatory formatting of meme input data that was generated from the updownstream.py script:
            - duplicate sample sequences must be dealt with to prevent MEME from biasing motifs from them
            - these duplicates arise in two main ways:
                1) Pulley samples:      genes with exactly the same 5' start share the upstream promoter    => delete one of the duplicates, the other gene gets nothing
                2) Tug'o'war samples:   a shared sample lies between opposite facing genes (+/-)            => half goes to one gene, half goes to the other

    Args:

    """
    # Need to ensure data is ready for MEME
    file_in  = open(filename_in)        # Read file in order to establish filters that will later format the data ready for MEME
    file_out = open(filename_out,'w')
    headers  = {}

    print '\tEstablishing filters...'
    while True:
        header  = file_in.readline()
        sequence= file_in.readline()
        if header == '':                    # break when finished
            break

    # tabs of the fi
        header_split    = header.split('\t')
        geneId          = header_split[0].replace('>','')
        utrCoord        = header_split[1].replace('UtrCoord:','').split(':')
        utrLength       = header_split[2].replace('UtrLength:','')
        utrType         = header_split[3].replace('UtrType:','').rstrip()
        chromosome      = utrCoord[0]

    # sub-properties of chromosome
        startEnd        = utrCoord[1].rpartition('-')   # rpartition ensures: (-800-20) -> '-800','-','20' 
        start           = startEnd[0]
        end             = startEnd[2]
        strand          = utrCoord[2]

    # dictionary storeage of all headers
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
    #   dealt with by making a set of unique utrcoords stripped of their strand property. Now upon encountering a match, forward samples and reverse samples take their half of the tug-o-war sequence
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
            if one_coord_per_geneId[':'.join(utrCoord)] == geneId:
                utrCoord_noStrand   = ':'.join(utrCoord[0:1+1])
            # TUG'O'WAR filter:     =>  Is this sample in a tug'o'war? 
                if utrCoord_noStrand in coords_shared:
                    coords_shared_geneIds.append(header.split('\t')[0].replace('>',''))
                    sequence = sequence[(len(sequence)/2):] # this works, since the rightermost end of the seq = gene, so we must sample from the end inwards
                # OKAY:             => then upon passing the filters, we can write the sample :)
                file_out.write(header+sequence+'\n')
    file_in.close()
    file_out.close()    

def memedata_dust( filename_in, cutoff=2 ): 
    """
    Notes: 
        DUSTMASKER: low complexity region masker, commandline program installed on vb-dev

    Args:
        cutoff = 5 <= the -level argument proportional to how strict the masking is, 20 is default for the native program, lower => more masking
    """
    file_dust_tmp_name  = filename_in.split('.')[0]+'_dusted_tmp'+'.'+filename_in.split('.')[1]
    file_dust_name      = filename_in.split('.')[0]+'_dusted'+'.'+filename_in.split('.')[1]
    print('\t\tMasking low complexity regions: DUST...')
    file_dust_tmp       = open(file_dust_tmp_name,'w')
    subprocess.call([ 'dust', filename_in, str(cutoff)],stdout=file_dust_tmp) # TODO: what value do I use for 'cut-off', here given as 10
    file_dust_tmp.close()
    file_dust_tmp   = open(file_dust_tmp_name)
    file_dust       = open(file_dust_name,'w')   # Concatenating consective lines of sequence...
    seq_united      = 'start'
    while True:
        line = file_dust_tmp.readline()
        if line == '':
            file_dust.write(seq_united)
            break
        if '>' in line:
            if seq_united == 'start':
                file_dust.write(line)
                seq_united = ''
            else:
                file_dust.write(seq_united+'\n'+line)
                seq_united = ''
        else:
            seq_united += line.rstrip()
    file_dust_tmp.close()
    file_dust.close()
    subprocess.call([ 'rm', file_dust_tmp_name]) # Deletes the temporary dust file

def memedata_simpleMask( filename_in ): 
    """
    Notes: 
        BOB's SIMPLE MASKER: low complexity region masker, masks di-nucleotide repeats if >=8bp or tri nucletides repeats if >=9bp long

    Args:

    """
    print('\t\tMasking low complexity regions: SIMPLE...')
    file_in  = open(filename_in)
    file_out = open('..'+filename_in.split('..')[1].split('.')[0]+'_simpleMasked'+'.'+filename_in.split('..')[1].split('.')[1],'w')
    def maskChars(matchobj):
         return len(matchobj.group(0))*'N'
    while True:
        line = file_in.readline()
        if line == "":  # break when finished
            break
        di_nucleotides = re.compile(r'(..)\1{4,}')      # Mask where DI NUCLEOTIDE repeats are 8bp long or more
        mask = di_nucleotides.sub(maskChars,line)
        tri_nucleotides = re.compile(r'(...)\1{3,}')    # Mask where TRI NUCLEOTIDE repeats are 9bp long or more
        mask = tri_nucleotides.sub(maskChars,mask)
        file_out.write(mask)
    file_in.close()
    file_out.close()

def allSpecies( speciesFile='./species_list.txt', dataPath_in='../data/sample_seqs/fasta/', dataPath_out = '../data/meme_data/in/', LCR_masking='simple' ): 
    """
    Notes:
        prepares the 2kb upstream data for MEME to start processing for every species

    Args:
        species_list = ['Anopheles gambiae', 'Aedes aegypti']    #list: species names to sample from
    """
    import speciesManage # species.generate_list() method creates a python list with names of species corresponding to MySQL EnsEMBL databases
    species_list = speciesManage.generate_list(speciesFile) # generates a list of species names corresponding to EnsEMBl MySQL db names

    for species in species_list:
        print(species)

        print(dataPath_in)

        filename_in = dataPath_in + species + '_upstream.fasta'
        filename_out= dataPath_out+ species + '_upstream_memeready_all.fasta'

        memedata_prep(filename_in,filename_out)

        if LCR_masking=='dust':
            memedata_dust(filename_out)
        elif LCR_masking=='simple':
            memedata_simpleMask(filename_out)

#-------------------------------------------------------------------------------------------------------------------------------------
# RUN
#-------------------------------------------------------------------------------------------------------------------------------------

""" NOTICE
     - dataPath_in  = .../masking/test/...
     - dataPath_out = .../meme_data/in/test...
"""

if __name__ == "__main__":
    import sys
    # Deal with MISSING ARGS:
    try: 
        speciesFile = sys.argv[1]
    except IndexError:
        speciesFile = './species_list.txt', # to alter the list of species to iterate through look in: scripts/species_list.txt
    try:
        dataPath_in = sys.argv[2]
    except IndexError:
        dataPath_in = '../data/sample_seqs/fasta/masking/',
    try:
        dataPath_out = sys.argv[3]
    except IndexError:
        dataPath_out = '../data/meme_data/in/'
    try:
        LCR_masking = sys.argv[4]
    except IndexError:
        LCR_masking = 'simple' # other options: 'dust', 'none'

    # RUN MAIN:
    allSpecies( speciesFile, dataPath_in, dataPath_out, LCR_masking )
    print('COMPLETE!')

# # PRODUCTION:
# allSpecies(    speciesFile  = './species_list.txt', # to alter the list of species to iterate through look in: scripts/species_list.txt
#                         dataPath_in  = 'home/ab108/0VB/2kb/data/sample_seqs/fasta/masking/',
#                         dataPath_out = 'home/ab108/0VB/2kb/data/meme_data/in/'
#                         )

# TG repeat story:
# allSpecies(    speciesFile  = './species_list.txt', # to alter the list of species to iterate through look in: scripts/species_list.txt
#                         dataPath_in  = 'home/ab108/0VB/2kb/data/sample_seqs/fasta/bob/',
#                         dataPath_out = 'home/ab108/0VB/2kb/data/meme_data/in/',
#                         dust=False
#                         )

# # DEBUGGING:
# allSpecies(    species_list, 
#                         dataPath_in  = 'home/ab108/0VB/2kb/data/sample_seqs/fasta/masking/test/',
#                         dataPath_out = 'home/ab108/0VB/2kb/data/meme_data/test/'
#                         )

