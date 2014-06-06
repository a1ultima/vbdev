
import re
import os 

#home = os.path.expanduser('~')
#os.chdir(home)

def get_repeat_geneIds(file_in):
    """
    Notes:

    Args:
        file_in     = '../data/meme_data/in/anopheles_gambiae_upstream_memeready_all_dusted.fasta'
        tg_repeats  = 3
    """
    file_in = open(file_in)                   # Open file for reading
    tg_query    = re.compile(r'(TG)\1{10,}')
    tg_samples  = []
    while True:
        header  = file_in.readline()
        seq     = file_in.readline()
        if header == "":                    # break when finished
            break
        if tg_query.search(seq):
            tg_samples.append(header.split('\t')[0].replace('>',''))
            #tg_locs.append(header.split('\t')[0].replace('>',''))
    file_in.close()
    return tg_samples

def gene_to_GO(file_in):
    """
    Notes:

    Args:
        file_in = open('0VB/2kb/data/meme_data/in/anopheles_gambiae_upstream_memeready_all_dusted.fasta')                   # Open file for reading
    """
    file_in = open(file_in)                   # Open file for reading
    tg_query    = re.compile(r'(TG)\1{3,}')
    tg_samples  = []
    while True:
        header  = file_in.readline()
        seq     = file_in.readline()
        if header == "":                    # break when finished
            break
        if tg_query.search(seq):
            tg_samples.append(header.split('\t')[0].replace('>',''))
    file_in.close()
    return tg_samples


# RUN:

if __name__ == "__main__":
    import sys
    # Deal with MISSING ARGS:
    try: 
        file_in = sys.argv[1]
    except IndexError:
        file_in = '../data/meme_data/in/anopheles_gambiae_upstream_memeready_all_dusted.fasta'

    # RUN MAIN:
    tg_samples = get_repeat_geneIds(file_in)
    print('COMPLETE!')



#file_in = '../data/tg_repeats/tg3_biomart_go.txt'



