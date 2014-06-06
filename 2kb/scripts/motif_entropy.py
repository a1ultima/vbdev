import numpy as np

def compute_entropy(motif_arr):
    
    H = -(  motif_arr[motif_arr > 0] * np.log2(motif_arr[motif_arr > 0])  ).sum(axis=0).mean() # WARN: not sure why .mean() is necessary... // SO: http://stackoverflow.com/questions/23480002/shannon-entropy-of-data-in-this-format-dna-motif // Formula from: http://en.wikipedia.org/wiki/Sequence_logo

    return H 

def fbp(fbp_path,format='decimal'):
    """ 
    Shannon entropy of a motif as float 

    Args:
        fbp_path    = '../data/stamp_data/out/dreme_100bp_e0.05/SWU_SSD/e005FBP.txt'
        format      = 'decimal' # also: 'count' // CTRL+F: "decimal" or "count"
    """
    fi      = open(fbp_path,'r')
    motif   = []
    while True:
        line = fi.readline().strip().lower()
        if line == "":
            break
        elif line.startswith('de'):
            pass
        elif line == 'xx': # when it's ended process the cached motif
            if motif:
                motif_arr = np.array(motif)
                if format == 'decimal':  
                    motif_H = compute_entropy(motif_arr)
                elif format == 'count':
                    motif_H = compute_entropy(motif_arr/motif_arr[0].sum())
            motif = []
        else: # append motif rows if the line is not a header, nor a delimiter, nor a break
            motif.append(map(float,line.split()[1:-1]))
    fi.close()
    return motif_H

if __name__ == "__main__":
    import sys
    try:
        fbp(sys.argv[1])
    except IndexError:
        fbp('../data/stamp_data/out/dreme_100bp_e0.05/SWU_SSD/e005FBP.txt')