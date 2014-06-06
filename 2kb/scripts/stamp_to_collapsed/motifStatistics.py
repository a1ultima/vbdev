
# MOTIF STATS
import os 
print(os.getcwd())
if not os.getcwd().endswith('scripts'):
    os.chdir('../')
    print os.getcwd()

import speciesManage
import numpy as np 
import matplotlib.pyplot as plt

save_d = 0.05
e = 0.05

filename_species= './species_list.txt'
species_list    = speciesManage.generate_list(filename_species) # generates a list of species names corresponding to EnsEMBl MySQL db name

# symbol_to_species dict    : to later translate motifs names into species that carry them
symbol_to_species = {}
for species in species_list:
    symbol = species.split('_')[0][0].upper()+species.split('_')[1][0:3].upper()
    symbol_to_species[symbol] = species

#------------------------------------------------------------------------------
# SPECIES -> MOTIFS 
#
# species_to_motifs         : to get #motifs per species
filename_dreme      = '../data/stamp_data/in/dreme_100bp_e'+str(e)+'/dreme.fasta'
file_dreme          = open(filename_dreme,'r') 
species_to_motifs   = {}
while True:
    line = file_dreme.readline()
    if line == "":               
        break
    elif line.startswith('DE '):               # Motif header
        motif  = line.split(' ')[1].rstrip()   # header (motif name)
        symbol = motif.split('_')[0]           # Species symbol
        # Collect a list of motifs per species
        if species_to_motifs.has_key(symbol_to_species[symbol]):
            species_to_motifs[symbol_to_species[symbol]].append(motif)
        else:
            species_to_motifs[symbol_to_species[symbol]] = [motif]
        while True:
            row = file_dreme.readline() # Row of pwm
            if ("XX" or "") in row:     # End of motif
                break
file_dreme.close()
# no. motifs per species
print 'motifs'
species_nmotifs_dist = []
for species in species_to_motifs.keys():
    nmotifs = len(species_to_motifs[species])
    species_nmotifs_dist.append(nmotifs)
    print('\t'+species + ': ' + str(len(species_to_motifs[species])))
species_nmotifs_dist = np.array(nmotifs)

#------------------------------------------------------------------------------
# SPECIES -> CLUSTERS
#
# no. clusters per species
#
clusters_path       = '../data/stamp_data/out/dreme_100bp_e'+str(e)+'/SWU_SSD/cluster_motifs_d'+str(save_d)+'/headers_for_fasta/'
clusters_files      = os.listdir(clusters_path)
species_to_clusters = {}
cluster_entropies   = []
for cluster_filename in clusters_files:
    cluster         = cluster_filename.replace('.txt','')
    cluster_file    = open(clusters_path+cluster_filename,'r')
# crawl through cluster file to read the motif names
    cluster_species = []
    while True:
        motif   = cluster_file.readline().rstrip()
        if motif == "":
            break
        symbol  = motif.split('_')[0]        
        species = symbol_to_species[symbol]
        cluster_species.append(species)
    cluster_file.close()

# append the cluster to species who were found in the cluster
    unique_cluster_species = set(cluster_species)    # get the unique list of species within cluster
    for species in unique_cluster_species:
        if species_to_clusters.has_key(species):
            species_to_clusters[species].append(cluster)
        else:
            species_to_clusters[species] = [cluster]
# no. clusters per species
print 'clusters'
species_nclusters_dist = []
for species in species_to_clusters.keys():
    nclusters = len(species_to_clusters[species])
    species_nclusters_dist.append(nclusters)
    print('\t'+species + ': ' + str(len(species_to_clusters[species])))
species_nclusters_dist = np.array(nclusters)

x   = range(0,len(species_list))
m   = species_nmotifs_dist
c   = species_nclusters_dist

# # diagram
# #plt.errorbar(x, m, xerr, yerr, capsize=0, ls='none', color='gray', elinewidth=1)

# # Text to annotate on plot
# m_labels = [str(int(round(i))) for i in m] # means 
# c_labels = [str(int(round(i))) for i in c] # medians

# # means
# for xpos, ypos, yname in zip(x, m, m_labels):
#     plt.annotate(yname, (xpos, ypos), xytext=(0, 0), va='bottom',textcoords='offset points')
# plt.annotate('Means (black)', (0.75, 50), xytext=(-10, 60), va='bottom',textcoords='offset points',arrowprops=dict(facecolor='black', shrink=0.05))

# # medians (red)
# for xpos, zpos, zname in zip(x, c, c_labels):
#     plt.annotate(zname, (xpos, zpos), xytext=(0, 0), va='bottom',textcoords='offset points',color='Red')
# plt.annotate('Medians (red)', (0.75, 15), xytext=(-10, -60), va='bottom',textcoords='offset points',color='Red',arrowprops=dict(facecolor='Red', shrink=0.05))

# # max 
# # for xpos, apos, aname in zip(x, a, adescrip):
# #     plt.annotate(aname, (xpos, apos), xytext=(0, 0), va='bottom',textcoords='offset points')
# # plt.annotate('Max (blue)', (0.75, 50), xytext=(-10, 60), va='bottom',textcoords='offset points',arrowprops=dict(facecolor='black', shrink=0.05))

# # quartiles
# # for xpos, upos, uname in zip(x, u, udescrip):
# #     plt.annotate(zname, (xpos, zpos), xytext=(0, 0), va='bottom',textcoords='offset points',color='Red')
# # plt.annotate('upper q', (0.75, 15), xytext=(-10, -60), va='bottom',textcoords='offset points',color='Red',arrowprops=dict(facecolor='Red', shrink=0.05))
# # for xpos, lpos, lname in zip(x, l, ldescrip):
# #     plt.annotate(zname, (xpos, zpos), xytext=(0, 0), va='bottom',textcoords='offset points',color='Red')
# # plt.annotate('lower q', (0.75, 15), xytext=(-10, -60), va='bottom',textcoords='offset points',color='Red',arrowprops=dict(facecolor='Red', shrink=0.05))

# # plotting
# plt.xlabel('Collapse clades with average distance to leaves < X (e.g. 0.2)')
# plt.ylabel('Mean size of clusters after collapse')
# plt.show()






