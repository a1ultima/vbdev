import speciesManage
import numpy as np 
import os 

suffix_in   = '_100bp/dreme.txt'
suffix_out  = '_100bp/dreme.fasta'

def allSpecies( speciesFile='./species_list.txt', dataPath_in='../data/meme_data/out/dreme_100bp/sampled_all_hpc/', method='STAMP' ):
    """
    Notes:
        - Converts motif formats from dreme -> e.g. STAMP (Transfac)
        - Unites all species motifs in one

    Args: 
        speciesFile = './species_list.txt'                                  # where is the species cfg file?
        dataPath_in = '../data/meme_data/out/dreme_100bp/sampled_all_hpc/'  # where are teh motifs?
        method      = 'STAMP'                                               # which format?
    """

    species_list                = speciesManage.generate_list(speciesFile) # generates a list of species names corresponding to EnsEMBl MySQL db names
    species_to_motifs_to_pwms   = {}

# Unite Motifs
    dataPath_out_united         = dataPath_in+'all_100bp'
    if not os.path.exists(dataPath_out_united):     # make a dir to store all species motifs in one transfac fasta
        os.makedirs(dataPath_out_united)
    filename_out_united = dataPath_in+'all'+suffix_out
    file_out_united     = open(filename_out_united,'w')

    for species in species_list:
    # READ MOTIF PWMs
        filename_in = dataPath_in+species+suffix_in
        file_in     = open(filename_in)
        species_to_motifs_to_pwms[species] = {}
        while True:
            line = file_in.readline()
            if line == "": 
                break
            if 'MOTIF ' in line:
                motif = line.split(' ')[1].rstrip()
                species_to_motifs_to_pwms[species][motif] = []
            if "letter-probability matrix" in line:
                stats   = line
                nsites  = int(line.split('=')[-2].replace(' ','').replace('E',''))
                pwm     = []
                while True:
                    line = file_in.readline()
                    if line == "\n":
                        species_to_motifs_to_pwms[species][motif] = {'pwm':np.array(pwm),'nsites':nsites,'stats':stats} #if we hit \n the files pwm is finished so we dict it
                        break
                    else:
                        pwm.append([float(i.rstrip()) for i in line.split(' ')]) #add the latest row of file's pwm to python pwm
        file_in.close()
    # WRITE FASTA PWMs
        filename_out= dataPath_in+species+suffix_out
        file_out    = open(filename_out,'w')
        for count,motif in enumerate(species_to_motifs_to_pwms[species].keys()):
            #header      = '>'+motif+'\n'
            pwm_counts  = np.round(species_to_motifs_to_pwms[species][motif]['pwm']*species_to_motifs_to_pwms[species][motif]['nsites']).astype(int) # convert PWM into strings
            pwm_check   = [ list(pwm_counts[i]) for i in range(0,len(pwm_counts)) ]
            pwm_str     = ''.join([' '.join([str(i) for i in row])+'\n' for row in pwm_check])

            if method=='STAMP':
                # ------------------------
                # Motif 2 position-specific probability matrix
                # ------------------------
                # letter-probability matrix: alength= 4 w= 6 nsites= 31
                # 0 31 0 0
                # 29 0 0 2
                # 0 30 0 1
                # 2 1 28 0
                # 0 3 0 28hhhhhh
                # 0 0 31 0
                sp      = species.split('_')[0][0].upper() + species.split('_')[1][0:3].upper()
                header  = 'Motif '+str(count+1)+' '+sp+'_'+motif+' position-specific probability matrix\n'
                stats   = species_to_motifs_to_pwms[species][motif]['stats']

                file_out.write(         '------------------------\n'+
                                        header.replace('>','')+ 
                                        '------------------------\n'+
                                        stats+
                                        pwm_str
                                        )
                file_out_united.write(  '------------------------\n'+
                                        header.replace('>','')+ 
                                        '------------------------\n'+
                                        stats+
                                        pwm_str
                                        )
            elif method=='MATLIGN':
                header      = '>'+species+'_'+motif+'\n'
                file_out.write(         header+
                                        pwm_str
                                        )
                file_out_united.write(  header+
                                        pwm_str
                                        )
            elif method=='TRANSFAC':
                # NA Mync
                # XX
                # DE Mync
                # XX
                # P0 A C G T
                # 01 0 31 0 0 C
                # 02 29 0 0 2 A
                # 03 0 30 0 1 C
                # 04 2 1 28 0 G
                # 05 0 3 0 28 T
                # 06 0 0 31 0 G
                # XX

                pwm_str = ''.join([str(j)+'\t'+str('\t'.join(i.split(' ')))+'\t'+motif[j]+'\n' for j,i in enumerate(pwm_str.split('\n')[:-1])]+['XX\n'])
                sp      = species.split('_')[0][0].upper() + species.split('_')[1][0:3].upper()
                print sp
                header  = 'DE '+sp+'_'+motif+'\n'
                file_out.write(         header+
                                        pwm_str
                                        )
                file_out_united.write(  header+
                                        pwm_str
                                        )
        file_out.close()
    file_out_united.close()
    return species_to_motifs_to_pwms

# RUN
if __name__ == "__main__":
    import sys
    try: 
        speciesFile = sys.argv[1]
    except IndexError:
        speciesFile = './species_list.txt'
    try: 
        dataPath_in = sys.argv[2]
    except IndexError:
        dataPath_in = dataPath_in = '../data/meme_data/out/dreme_100bp/sampled_all_hpc/evalue_0.05/'
    try:
        method = sys.argv[3]
    except IndexError:
        method = 'TRANSFAC'
    species_to_motifs_to_pwms = allSpecies( speciesFile, dataPath_in, method)
    print('COMPLETE!')


#################
#   Example     #
#################

# From iPython:
# import motifPmw_toFasta
# cool = motifPmw_toFasta.allSpecies( speciesFile='./species_list.txt', dataPath_in='../data/meme_data/out/dreme_100bp/sampled_all_hpc/evalue_0.05/', method='TRANSFAC' )

