import speciesManage

speciesFile = './species_list.txt'
dataPath_in = '../data/meme_data/out/dreme_100bp/sampled_all_hpc/'
suffix_in   = '_100bp/dreme.txt'
species_to_motifs = {}
species_list = speciesManage.generate_list(speciesFile)

# READ DREME => species:motifs
for species in species_list:
    print(species)
    species_to_motifs[species]= []
    filename_in = dataPath_in+species+suffix_in
    # Read unique motifs
    try:
        file_in = open(filename_in,'r')
        count = 0
        while True:
            line = file_in.readline()
            if line == "":
                break
            if 'MOTIF' in line:
                count += 1
                motif = line.split(' ')[1].rstrip()   # split line where '\t' occurs
                species_to_motifs[species].append(motif)
    except IOError:
        file_in.close()
        print('skipping: '+species)
        continue
    file_in.close()

# READ TOMTOM => species:species:matchedMotifs
dataPath_in = '../data/meme_data/out/tomtom_100bp/'
suffix_in   = '_100bp/tomtom.txt'
speciesPair_to_motifMatches = {}
for species_query in species_list:
    print(species_query)
    speciesPair_to_motifMatches[species_query] = {}
    for species_database in species_list:
        print('\t'+species_database)
        speciesPair_to_motifMatches[species_query][species_database] = {'list':[],'set':[]}
        filename_in  = dataPath_in+species_query+'_vs_'+species_database+suffix_in
        file_in      = open(filename_in,'r')
        motif_queries= []
        while True:
            line = file_in.readline()
            if line == "":
                break
            if line.startswith('#'): # skip header
                continue
            line_split = line.split('\t')
            motif_query = line_split[0]     # the query motif that is matched
            motif_queries.append(motif_query)
            speciesPair_to_motifMatches[species_query][species_database]['list'] = motif_queries
            speciesPair_to_motifMatches[species_query][species_database]['set']  = set(motif_queries)
        file_in.close()

#---------------------------------------------------------------------------------------------

#nMotifs_distribution = sorted([(len(species_to_motifs[species]),species) for species in species_to_motifs.keys()])
nMotifs_distribution = sorted([len(species_to_motifs[species]) for species in species_to_motifs.keys()])

import pickle
pickle.dump(nMotifs_distribution,open('nMotifs_distribution.p','wb'))
pickle.dump(speciesPair_to_motifMatches,open('speciesPair_to_motifMatches.p','wb'))
pickle.dump(species_to_motifs,open('species_to_motifs.p','wb'))

# # a bar plot with errorbars
# import numpy as np
# import matplotlib.pyplot as plt

# N = 5
# menMeans = nMotifs_distribution
# menStd =   [1]*len(nMotifs_distribution)

# ind = np.arange(N)  # the x locations for the groups
# width = 0.35       # the width of the bars

# fig, ax = plt.subplots()
# rects1 = ax.bar(ind, menMeans, width, color='r', yerr=menStd)

# womenMeans = nMotifs_distribution
# womenStd =   [1]*len(nMotifs_distribution)
# rects2 = ax.bar(ind+width, womenMeans, width, color='y', yerr=womenStd)

# # add some
# ax.set_ylabel('Scores')
# ax.set_title('Scores by group and gender')
# ax.set_xticks(ind+width)
# ax.set_xticklabels( ('G1', 'G2', 'G3', 'G4', 'G5') )

# ax.legend( (rects1[0], rects2[0]), ('Men', 'Women') )

# def autolabel(rects):
#     # attach some text labels
#     for rect in rects:
#         height = rect.get_height()
#         ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
#                 ha='center', va='bottom')

# autolabel(rects1)
# autolabel(rects2)

# plt.show()