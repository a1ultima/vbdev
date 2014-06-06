import subprocess 
import os 

def meme_bfileGenerator( file_in, file_out, bfileGenerator_path, masking_char = 'n', order=3 ):
    """
    Notes:
        Generates MEME bacground markov models of order 2 which is used by MEME in the -bfile arg 

    Dependencies:
        Executed:   updownstream.py --> meme_dataPrepper.py --> meme_bfileGenerator.py
        Scripts:    Bob's perl bfileGenerator.pl script 

    Args: 
    """
    print('\tgenerating bfile...')
    home = os.path.expanduser('~')  # set path to home for absolute pathing later
    os.chdir(home)
# BOB's BFILE GENERATOR .pl script
    bfileGenerator_cmd = [  'perl', 
                            bfileGenerator_path+'bfileGenerator.pl', 
                            '-order', str(order), 
                            file_in
                            ]
    file_out_tmp_name = file_out.replace('.bg2','_tmp.bg2')
    subprocess.call( bfileGenerator_cmd,stdout=open(file_out_tmp_name,'w') ) # generte bfile
# PRUNE OUT MASKING 'N's    <= meme prohibits 'N's in the -bfile, so we must prune em out
    file_in_tmp = open(file_out_tmp_name,'r')   # Open file for reading
    file_out    = open(file_out,'w')            # Open file for writing
    print '\tpruning illegal characters...'
    while True:
        line = file_in_tmp.readline()
        if line == "":                    # break when finished
            break
        if not masking_char in line:
            file_out.write(line)
        else:
            continue
    file_in_tmp.close()
    file_out.close()
    subprocess.call([ 'rm', file_out_tmp_name]) # Deletes the temporary illegal bfile file

def allSpecies( inputData_path, bfileGenerator_path, speciesFile='./species_list.txt', masking_char = 'n', order=3 ):
    """
    Notes:
        Iterates meme_bfileGenerator for all species listed

    Dependencies:
        Executed:   see meme_bfileGenerator.py 
        Scripts:    see meme_bfileGenerator.py

    Args: 
        inputData_path = '0VB/2kb/data/meme_data/in/'      # path of directory with allSpecies memeready fasta files
    """
    print('Generating bfiles for all species...')
    import speciesManage # species.generate_list() method creates a python list with names of species corresponding to MySQL EnsEMBL databases
    species_list = speciesManage.generate_list(speciesFile) # generates a list of species names corresponding to EnsEMBl MySQL db names
    # species_list = [    'Aedes aegypti',
    #                 'Anopheles albimanus',
    #                 'Anopheles arabiensis',
    #                 'Anopheles atroparvus',
    #                 'Anopheles christyi',
    #                 'Anopheles culicifacies',
    #                 'Anopheles darlingi',
    #                 'Anopheles dirus',
    #                 'Anopheles epiroticus',
    #                 'Anopheles farauti',
    #                 'Anopheles funestus',
    #                 'Anopheles gambiae',
    #                 'Anopheles maculatus',
    #                 'Anopheles melas',
    #                 'Anopheles merus',
    #                 'Anopheles minimus',
    #                 'Anopheles quadriannulatus',
    #                 'Anopheles sinensis',
    #                 'Anopheles stephI',           # this is the renamed stephensiI db
    #                 'Anopheles stephensi'
    #                 ]
    for species in species_list:
        print(species)
        #file_in     = inputData_path+species.lower().replace(' ','_')+'_upstream_memeready_all_dusted.fasta'   # ANDY 25_02
        file_in     = inputData_path+species.lower().replace(' ','_')+'_upstream_memeready_all.fasta'           # ANDY 25_02
        file_out    = inputData_path+species.lower().replace(' ','_')+'.bg2'
        meme_bfileGenerator( file_in, file_out, bfileGenerator_path, masking_char, order ) 

#-------------------------------------------------------------------------------------------------------------------------------------
# RUN
#-------------------------------------------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    import sys
# Deal with MISSING ARGS:
    try: 
        inputData_path = sys.argv[1]
    except IndexError:
        inputData_path = '0VB/2kb/data/meme_data/in/'  # where to find .fasta files for MEME to work on
    try:
        bfileGenerator_path = sys.argv[2]
    except IndexError:
        bfileGenerator_path = '0VB/2kb/scripts/'            # where to find Bob's perl bfile script
    try:
        masking_char = sys.argv[3]
    except IndexError:
        masking_char = 'n'
    try:
        order = sys.argv[4]
    except IndexError:
        order = 3
# RUN:
    allSpecies( inputData_path, bfileGenerator_path, masking_char = 'n', order=3, speciesFile='./species_list.txt' )
    print('COMPLETE!')