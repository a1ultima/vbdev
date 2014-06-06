
# MOTIF STATS
import os 
print(os.getcwd())
if not os.getcwd().endswith('scripts'):
    os.chdir('../')
    print os.getcwd()

import speciesManage
import motif_entropy

import numpy as np 
#import matplotlib.pyplot as plt

save_d  = 0.05
e       = 0.05

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

clusters_fbp_path       = '../data/stamp_data/out/dreme_100bp_e'+str(e)+'/SWU_SSD/cluster_motifs_d'+str(save_d)+'/family_binding_matrices/'
clusters_tree_path      = '../data/stamp_data/out/dreme_100bp_e'+str(e)+'/SWU_SSD/cluster_motifs_d'+str(save_d)+'/family_binding_matrices/'
clusters_jaspar_path    = '../data/stamp_data/out/dreme_100bp_e'+str(e)+'/SWU_SSD/cluster_motifs_d'+str(save_d)+'/family_binding_matrices/'
clusters_transfac_path  = '../data/stamp_data/out/dreme_100bp_e'+str(e)+'/SWU_SSD/cluster_motifs_d'+str(save_d)+'/family_binding_matrices/'
clusters_path           = '../data/stamp_data/out/dreme_100bp_e'+str(e)+'/SWU_SSD/cluster_motifs_d'+str(save_d)+'/headers_for_fasta/'

clusters_files      = os.listdir(clusters_path)
species_to_clusters = {}
cluster_to_stats    = {}
cluster_entropies   = []

for cluster_filename in clusters_files:
    
    cluster = cluster_filename.replace('.txt','')    
    
    cluster_to_stats[cluster] = {
        'species':{
            'unique':{
                'n_unique':0,
                'n_paralogue':0,
                'list':[]
                },
            'phy_profile':{
                'binary':{},
                'count':{}
                }
            },
        'motif':{
            'n':0,
            'list':[]
            },
        'cluster':{
            'path':{
                'fbp':str(),
                'tree':str(),
                'jaspar':str(),
                'transfac':str()},
            'data':{
                'fbp':0,
                'tree':0,
                'jaspar':0,
                'transfac':0},
            'H':float()
            }
        }

    for species in species_list:
        cluster_to_stats[cluster]['species']['phy_profile']['binary'][species] = 0
        cluster_to_stats[cluster]['species']['phy_profile']['count'][species] = 0

# Initiate filenames

    # Filenames  
    cluster_tree    = clusters_tree_path+cluster+'.tree'
    cluster_fbp     = clusters_fbp_path+cluster+'FBP.txt'
    cluster_jaspar  = clusters_jaspar_path+cluster+'_match_paris.txt'
    cluster_transfac= clusters_transfac_path+cluster+'_matched.transfac'

    cluster_to_stats[cluster]['cluster']['path']['tree']    = cluster_tree
    cluster_to_stats[cluster]['cluster']['path']['fbp']     = cluster_fbp
    cluster_to_stats[cluster]['cluster']['path']['jaspar']  = cluster_jaspar
    cluster_to_stats[cluster]['cluster']['path']['transfac']= cluster_transfac

    
# CLUSTER HEADERS           : crawl through cluster file to read the motif names
    cluster_file    = open(clusters_path+cluster_filename,'r')
    cluster_species = []
    while True:
        motif   = cluster_file.readline().rstrip()
        if motif == "":
            break
        symbol  = motif.split('_')[0]        
        species = symbol_to_species[symbol]
        cluster_species.append(species)
        cluster_to_stats[cluster]['motif']['n'] = 1 + cluster_to_stats[cluster]['motif']['n']
        cluster_to_stats[cluster]['motif']['list'].append(motif)
    cluster_file.close()

# FAMILY BINDING PROFILES   : crawl through cluster_fbp file to calculate the entropy

    cluster_fbp_entropy                         = motif_entropy.fbp(cluster_to_stats[cluster]['cluster']['path']['fbp'])
    cluster_to_stats[cluster]['cluster']['H']   = cluster_fbp_entropy # calculate entropy for the motif

    #cluster_tree.close()
    #cluster_fbp.close()
    #cluster_jaspar.close()
    #cluster_transfac.close()

# append the cluster to species who were found in the cluster
    for species in cluster_species:
        cluster_to_stats[cluster]['species']['phy_profile']['binary'][species] = 1
        cluster_to_stats[cluster]['species']['phy_profile']['count'][species] =  1 + cluster_to_stats[cluster]['species']['phy_profile']['count'][species]
    unique_cluster_species = set(cluster_species)    # get the unique list of species within cluster

    # cluster_to_species
    cluster_to_stats[cluster]['species']['unique']['n_unique'] = len(unique_cluster_species)
    cluster_to_stats[cluster]['species']['unique']['n_paralogue'] = len(cluster_species) - len(unique_cluster_species)
    cluster_to_stats[cluster]['species']['unique']['list'] = unique_cluster_species

    for species in unique_cluster_species:
        # species_to_cluster
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
#c   = species_nclusters_dist

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TODO: read the stamp data instead of just path names
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

clusters_path       = '../data/stamp_data/out/dreme_100bp_e'+str(e)+'/SWU_SSD/cluster_motifs_d'+str(save_d)+'/headers_for_fasta/'

event_to_cluster = {'e0.05_d0.05':{},
                    'e0.5_d0.05':{}
                    }

# clusters w/ species >= X
for s in range(1,len(species_list)): # per species number, S
    # find clusters with nSpecies >= s
    clusters_s = [c for c in cluster_to_stats.keys() if cluster_to_stats[c]['species']['unique']['n_unique'] >= s]
    
    # append to event_to_cluster dict
    if not event_to_cluster['e0.05_d0.05'].has_key('nSpeciesOrLess'):
        event_to_cluster['e0.05_d0.05'] = {'nSpeciesOrLess':{s:{'list':clusters_s,'nClusters':len(clusters_s),'H':{'distribution':[],'mean':float()}}}}
    else:
        event_to_cluster['e0.05_d0.05']['nSpeciesOrLess'][s] = {'list':     clusters_s, # same as above but once the 'nSpeciesOrLess' key is available
                                                                'nClusters':len(clusters_s),
                                                                'H':   float()}

# all clusters      : I put this after the species >=X since the "if not" check above can delete previously existing keys, I should use "update"
clusters_all = cluster_to_stats.keys()
event_to_cluster['e0.05_d0.05']['all'] = {
    'list':         clusters_all,
    'nClusters':    len(clusters_all),
    'H':       float()
}


# clusters w/ drosophilla (dipteran)
clusters_dmel = [c for c in cluster_to_stats.keys() if 'droso' in ''.join(cluster_to_stats[c]['species']['unique']['list'])]

# if not event_to_cluster['e0.05_d0.05'].has_key('inc_DMEL'):
#     event_to_cluster['e0.05_d0.05'] = {'inc_DMEL':{ 'list':clusters_dmel,
#                                                 'nClusters':len(clusters_dmel),
#                                                 'H':float()}}    

event_to_cluster['e0.05_d0.05']['inc_DMEL'] = { 
    'list':     clusters_dmel,
    'nClusters':len(clusters_dmel),
    'H':   float()
    }

#cluster_to_stats[cluster]['cluster']['H']

for cluster_set in event_to_cluster['e0.05_d0.05'].keys():
    
    clusters            = event_to_cluster['e0.05_d0.05'][cluster_set]['list']
    cluster_entropies   = np.array([])

    for cluster in clusters:
        np.append(cluster_entropies,cluster_to_stats[cluster]['cluster']['H'])








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


