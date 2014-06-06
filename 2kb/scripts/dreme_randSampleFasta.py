import sys 
import os 
import random
import re
import speciesManage # species.generate_list() method creates a python list with names of species corresponding to MySQL EnsEMBL databases

def allSpecies(speciesFile, dataPath_in, dataPath_out, n_seqs = 10000, len_seq = 100, population = True):
    """
    Notes: 

    Args:
        speciesFile = './species_list.txt'                  # location of species config file
        dataPath_in = '../data/meme_data/in/'               # input directiory
        dataPath_out='../data/meme_data/in/random_dreme/'   # output directory
        n_replicates=5                                      # number of replicates per species
        n_seqs      =10000                                  # 10k sampled seqs
        len_seq     =100                                    # 100bp long sequence samples
    """
    def rand_parts(seq, n, l):
        indices = xrange(len(seq) - (l - 1) * n)
        result = []
        offset = 0
        for i in sorted(random.sample(indices, n)):
            i += offset
            result.append(seq[i:i+l])
            offset += l - 1
        return result
    def chunks(l, n):
        """ Yield successive n-sized chunks from l. """
        for i in xrange(0, len(l), n):
            yield l[i:i+n]

        
    sys.path.append(os.getcwd())
    suffix_in       = '_upstream_memeready_all_simpleMasked.fasta'  # common filename component between al sp.
    suffix_out      = suffix_in.split('.')[0].replace('memeready','dremeready')+'_random' +'.'+suffix_in.split('.')[1]
    species_list    = speciesManage.generate_list(speciesFile) # generates a list of species names corresponding to EnsEMBl MySQL db names
    for species in species_list:
        print(species)
    # IMPORT FASTA                      : 
        filename_in = dataPath_in+species+suffix_in
        filename_out= dataPath_out+species+suffix_out
        file_in     = open(filename_in)
        headers     = {} # dict for geneId:header
        seqs        = {} # geneId:seq
        print '\tImporting FASTA...'
        while True:
            header  = file_in.readline()
            seq     = file_in.readline()
            if header == '':
                break
            header_split    = header.split('\t')
            geneId          = header_split[0].replace('>','')
            mask = re.compile(r'(N)\1{2,}')  # collapse consecutive masking chars, NNN --> N
            seq_collaspedMasks = mask.sub('N',seq).rstrip()+'N'  
            #                                                ^ add X to delimit each sequence later when they are concated
            headers[geneId] = header 
            seqs[geneId]    = seq_collaspedMasks
        file_in.close()
    # META SEQUENCE                     : concatenate all sequences for later random sampling
        print('\tconcatenating sequences to form meta-sequence...')
        meta_seq = ''.join([seqs[i] for i in seqs.keys()])
        meta_index = range(0,len(meta_seq))

    # LOCATION:GENEID                   : map each bp location of the meta_sequence to geneId it belongs to, ~14mil long meta_seq becomes its own key:geneId
        print('\tmapping meta-sequence to geneIds...')
        loc2GeneId  = {} 
        lens        = [len(seqs[i]) for i in seqs.keys()] # lengths of each seq 
        start       = 0 
        for k,geneId in enumerate(seqs.keys()):
            start   = start
            end     = start + lens[k]
            seqRange= range(start,end)
            for loc in seqRange:
                loc2GeneId[loc] = geneId
            start   = end
    # RANDOM NON_OVERLAPPING SAMPLING   : sample N random subranges along the meta_seq, representing non-overlapping sub_seqs length L
        print('\tgenerating random sub-sequence samples...')
        if population == False:
            randRanges = rand_parts(meta_index,n_seqs,len_seq)
        elif population == True:
            randRanges = chunks(meta_index,len_seq)
        randSamples = []
        for randRange in randRanges:
            geneSet = list(set([loc2GeneId[loc] for loc in randRange])) # collect geneIds whose sequences lie in the randRange 
            randSamples.append({'geneId':geneSet, # for the given rangeRange take note of the geneIds within
                                'seq':meta_seq[randRange[0]:randRange[-1]] # store the sub_seq corresponding to the randRange
                                })
    # GENERATE FASTA
        print('\twriting sampled sequences to FASTA...')

        if not os.path.exists(dataPath_out): # generate cluster fasta dirs
            print('making directory: '+dataPath_out)
            os.makedirs(dataPath_out)

        file_out = open(filename_out,'w')
        for sample in randSamples:
            headers = sample['geneId'] 
            seq = sample['seq']
            file_out.write('>'+'\t'.join(headers)+'\n')
            file_out.write(seq+'\n')
        file_out.close()

#----------------------------------------------------------------------

if __name__ == "__main__":
    try: 
        speciesFile = sys.argv[1]
    except IndexError:
        speciesFile = './species_list.txt'
    try: 
        dataPath_in = sys.argv[2]
    except IndexError:
        dataPath_in = '../data/meme_data/in/'
    try: 
        dataPath_out = sys.argv[3]
    except IndexError:
        dataPath_out = '../data/meme_data/in/random_dreme/'
    try:
        n_seqs = int(sys.argv[4])
    except IndexError:
        n_seqs = 10000
    try: 
        len_seq = int(sys.argv[5])
    except IndexError:
        len_seq = 100
    try: 
        population = bool(sys.argv[6])
    except IndexError:
        population = True
    allSpecies( speciesFile, dataPath_in, dataPath_out, n_seqs, len_seq, population )
    print('COMPLETE!')