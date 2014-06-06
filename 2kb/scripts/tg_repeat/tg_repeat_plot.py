

import re
import os 

home = os.path.expanduser('~')
os.chdir(home)

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
                'Anopheles sinensis',
                'Anopheles stephI',           # this is the renamed stephensiI db
                'Anopheles stephensi'
                ]

file_prefix = '0VB/2kb/data/meme_data/in/'

bfile_dict = {}

for species in species_list:

    species = species.lower().replace(' ','_')
    bfile_dict[species] = {}

    file_in = file_prefix+species+'.bg2'
    file_in = open(file_in)                   # Open file for reading

    bfile_dist = {}

    while True:
        line  = file_in.readline()
        if line == "":                    # break when finished
            break
        elif '#' in line:
            continue
        
        line_split = line.split(' ')

        motif   = line_split[0]
        freq    = line_split[-1]
        
        bfile_dict[species][motif] = float(freq)

    file_in.close()


freq_list = []
splist = []

for species in species_list:
    species = species.lower().replace(' ','_')
    splist.append(species)
    freq_list.append(bfile_dict[species]['tgtgtg'])

#import numpy as np 
# import pylab as pl

# fig = pl.figure()
# ax = pl.subplot(111)
# ax.bar(bfile_dict.keys(), freq_list, width=100)
# ax.xaxis_date()




