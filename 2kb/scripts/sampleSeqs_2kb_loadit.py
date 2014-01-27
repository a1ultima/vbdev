# IMPORTS
import pickle
import os
from cogent import LoadSeqs, DNA

#import matplotlib.pyplot as plt    # TODO: server backend needs install: 'ImportError: Gtk* backend requires pygtk to be installed.'


""" Convert sequences into a cogent SequenceCollection opbject
"""
print('gathering genes for fasta...')

pathin              = '../data/'     # path of the batched sample data
sample_filenames    = os.listdir(pathin)      # list all the sample batch file names
samples             = {}
gene_lengths        = []

for filename in sample_filenames:
    sample_data = pickle.load(open(pathin+filename,'rb')) # Load the sequences from pickled batches 
    geneIds = sample_data.keys() # gene ids 
    for gene in geneIds:        # dictionary for sequences
        if not samples.has_key(gene):
            if sample_data[gene]['truncated']:
                print('\t'+gene)
                gene_length = len(sample_data[gene]['truncated'])
                gene_lengths.append(gene_length)
                samples[gene+'_'+str(gene_length)+'bp'] = sample_data[gene]['truncated']
            else:
                print('\t'+gene+' (no sequence)')
        else:
            print('\t\tduplciate sample: '+gene)
            break

samples_loaded  = LoadSeqs(data=samples,moltype=DNA,aligned=False)
samples_fasta   = sampels_loaded.toFasta()

# import pprint
# pprint.pprint(sorted(gene_lengths))
#plt.hist(gene_lengths)
#plt.show()         # TODO: server backend needs install: 'ImportError: Gtk* backend requires pygtk to be installed.'



fileout_fasta = open(pathin+'all.fasta','wb')
fileout_fasta.write(samples_fasta)
fileout_fasta.close()