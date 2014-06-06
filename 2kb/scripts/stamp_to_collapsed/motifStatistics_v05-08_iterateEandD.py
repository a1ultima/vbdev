
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


event_to_clusterClass   = {}                        # events as primary key
# parameterisation_events = [(0.05,0.05),(0.5,0.5)]   # parameterisations to draw data from
parameterisation_events = [(0.05,0.05)]   # parameterisations to draw data from

for e,save_d in parameterisation_events:

    print('Parameters: '+str(e)+' '+str(save_d))

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
    clusters_fasta_path     = '../data/stamp_data/out/dreme_100bp_e'+str(e)+'/SWU_SSD/cluster_motifs_d'+str(save_d)+'/fasta/'
    clusters_path           = '../data/stamp_data/out/dreme_100bp_e'+str(e)+'/SWU_SSD/cluster_motifs_d'+str(save_d)+'/headers_for_fasta/'

    clusters_files          = os.listdir(clusters_path)
    species_to_clusters     = {}
    cluster_to_stats        = {}
    cluster_entropies       = []

    for count,cluster_filename in enumerate(clusters_files):
        
        cluster = cluster_filename.replace('.txt','')    
        
        #print(count,cluster)

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
            cluster_to_stats[cluster]['species']['phy_profile']['binary'][species]  = 0
            cluster_to_stats[cluster]['species']['phy_profile']['count'][species]   = 0

    # Initiate filenames
        # Filenames  
        cluster_tree    = clusters_tree_path+cluster+'.tree'
        cluster_jaspar  = clusters_jaspar_path+cluster+'_match_paris.txt'
        cluster_transfac= clusters_transfac_path+cluster+'_matched.transfac'
        cluster_fbp     = clusters_fbp_path+cluster+'FBP.txt'
        
        if not cluster_fbp in ''.join(os.listdir(clusters_path)):  # check if a FBP is available for this cluster (i.e. if nMotifs in the cluster is = 1 then there will not be a FBP due to stamp needing at least two motifs), if not available then use the .fasta instead
            cluster_fbp = clusters_fasta_path+cluster+'.fasta'

        cluster_to_stats[cluster]['cluster']['path']['tree']    = cluster_tree
        cluster_to_stats[cluster]['cluster']['path']['fbp']     = cluster_fbp

        cluster_to_stats[cluster]['cluster']['path']['jaspar']  = cluster_jaspar
        cluster_to_stats[cluster]['cluster']['path']['transfac']= cluster_transfac

    # CLUSTER HEADERS           : crawl through cluster file to read the motif names
        cluster_file    = open(clusters_path+cluster_filename,'r')
        cluster_species = []
        while True:
            motif = cluster_file.readline().rstrip()
            if motif == "":
                break
            symbol  = motif.split('_')[0]        
            species = symbol_to_species[symbol]
            cluster_species.append(species)
            cluster_to_stats[cluster]['motif']['n'] = 1 + cluster_to_stats[cluster]['motif']['n']
            cluster_to_stats[cluster]['motif']['list'].append(motif)
        cluster_file.close()

    # FAMILY BINDING PROFILES   : crawl through cluster_fbp file to calculate the entropy

        # try: # by default we apply the decimal format (works for clusters >= 2 motifs), if that errors we use the counts format (clusters = 1 motif)
        #     cluster_fbp_entropy                         = motif_entropy.fbp(cluster_to_stats[cluster]['cluster']['path']['fbp'],format='decimal')
        #     cluster_to_stats[cluster]['cluster']['H']   = cluster_fbp_entropy # calculate entropy for the motif
        # except IOError:
        #     cluster_fbp_entropy                         = motif_entropy.fbp(cluster_to_stats[cluster]['cluster']['path']['fbp'],format='count')
        #     cluster_to_stats[cluster]['cluster']['H']   = cluster_fbp_entropy

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

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 
    # TODO: read the stamp data instead of just path names # 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

    clusters_path       = '../data/stamp_data/out/dreme_100bp_e'+str(e)+'/SWU_SSD/cluster_motifs_d'+str(save_d)+'/headers_for_fasta/'

    # EVENT-TO-CLUSTER      :  an event is the data given a dreme evalue cut-off and a distance cut-off for tree-collapsing // common to each event are classes of interesting clusters, e.g. those with >= 3 species, or those with at least one drosophila motif, etc.// common to each cluster class are three pieces of key-value data: (i) entropy dict (mean, distribution, etc), (ii) list of cluster names in that class, (iii) and the number fo clusters. Due to common cluster classes per event and common cluster key-value data per cluster class we can construct the dictionary using "template dicts" for the cluster classes and for the key-value data.
    #  

    #------------------------------------------#
    #  nMotifs >= S, nMotifs == dipterans      #
    #------------------------------------------#

    event_to_clusterClass['e'+str(e)+'_d'+str(save_d)] = {}

    # events mapped to common clusterClasses
    event_to_clusterClass_template = { 
        'all':{},
        'nSpeciesOrMore':{},
        'clades':{
            'inc_DMEL':{},
            'inc_gambiaeComplex':{}
            }
        }
    # cluster classes mapped to common data
    clusterClass_to_stats_template = {
        'H':{
            'distribution':[],
            'mean':float()
            },
        'list':[],
        'nClusters':int()
        } 

    # #------------------------------------------#

    # Build 
    for event in event_to_clusterClass.keys():
        
        print event

        event_to_clusterClass[event] = event_to_clusterClass_template.copy() # cluster classes

        # all               :   summary stats for all the clusters
        clusters_all                                                = cluster_to_stats.keys()
        clusters_all_entropy                                        = np.array([cluster_to_stats[c]['cluster']['H'] for c in clusters_all]) # entropies for each cluster
        event_to_clusterClass[event]['all']                         = clusterClass_to_stats_template.copy()
        event_to_clusterClass[event]['all']['list']                 = clusters_all
        event_to_clusterClass[event]['all']['nClusters']            = len(clusters_all)
        event_to_clusterClass[event]['all']['H']['distribution']    = clusters_all_entropy
        event_to_clusterClass[event]['all']['H']['mean']            = clusters_all_entropy.mean()

        # nSpeciesOrMore    :   clusters w/ species >= X
        for s in range(1,len(species_list)+1): # per species number, S
            
            clusters_s                                                              = [c1 for c1 in cluster_to_stats.keys() if cluster_to_stats[c1]['species']['unique']['n_unique'] >= s] # find clusters with nSpecies >= s
            clusters_s_entropy                                                      = np.array([cluster_to_stats[c]['cluster']['H'] for c in clusters_s]) # entropies for each cluster

            event_to_clusterClass[event]['nSpeciesOrMore'][s]                       = clusterClass_to_stats_template.copy()
            event_to_clusterClass[event]['nSpeciesOrMore'][s]['list']               = clusters_s
            event_to_clusterClass[event]['nSpeciesOrMore'][s]['nClusters']          = len(clusters_s)

            event_to_clusterClass[event]['nSpeciesOrMore'][s]['H']['distribution']  = clusters_s_entropy
            event_to_clusterClass[event]['nSpeciesOrMore'][s]['H']['mean']          = clusters_s_entropy.mean()

        print(len(event_to_clusterClass[event]['all']['list']))

        # clades            :   one level deeper than the above, under clades we have various clade related conditions to define the classes the template must be appended to each
        for clade in event_to_clusterClass_template['clades'].keys():
            event_to_clusterClass[event]['clades'][clade] = clusterClass_to_stats_template.copy() # append common data template to each clade

        # inc_DMEL          :   clusters including at least one drosophila
        clusters_inc_DMEL                                                       = [c for c in cluster_to_stats.keys() if 'droso' in ''.join(cluster_to_stats[c]['species']['unique']['list'])] # find clusters that include drosophila motif/s

        clusters_inc_DMEL_entropy                                               = np.array([cluster_to_stats[c]['cluster']['H'] for c in clusters_inc_DMEL]) # entropies for each cluster

        event_to_clusterClass[event]['clades']['inc_DMEL']['list']              = clusters_inc_DMEL
        event_to_clusterClass[event]['clades']['inc_DMEL']['nClusters']         = len(clusters_inc_DMEL)
        event_to_clusterClass[event]['clades']['inc_DMEL']['H']['distribution'] = clusters_inc_DMEL_entropy
        event_to_clusterClass[event]['clades']['inc_DMEL']['H']['mean']         = clusters_inc_DMEL_entropy.mean()

        #clusters_inc_gambiaeComplex = 

    # print(len(event_to_clusterClass['e0.05_d0.05']['all']['list']))


    # event_to_clusterClass['e0.05_d0.05']['nSpeciesOrMore'] = {}

    # for s in range(1,len(species_list)+1): # per species number, S
        
    #     clusters_s                                                              = [c for c in cluster_to_stats.keys() if cluster_to_stats[c]['species']['unique']['n_unique'] >= s] # find clusters with nSpecies >= s
    #     clusters_s_entropy                                                      = np.array([cluster_to_stats[c]['cluster']['H'] for c in clusters_s]) # entropies for each cluster
    #     event_to_clusterClass['e0.05_d0.05']['nSpeciesOrMore'][s]                       = clusterClass_to_stats_template
    #     event_to_clusterClass['e0.05_d0.05']['nSpeciesOrMore'][s]['list']               = clusters_s
    #     event_to_clusterClass['e0.05_d0.05']['nSpeciesOrMore'][s]['nClusters']          = len(clusters_s)
    #     event_to_clusterClass['e0.05_d0.05']['nSpeciesOrMore'][s]['H']['distribution']  = clusters_s_entropy
    #     event_to_clusterClass['e0.05_d0.05']['nSpeciesOrMore'][s]['H']['mean']          = clusters_s_entropy.mean()

    # species_diversity_cluster_distribution = sorted([cluster_to_stats[c]['species']['unique']['n_unique'] for c in cluster_to_stats.keys()])

    # s = [event_to_clusterClass['e0.05_d0.05']['nSpeciesOrMore'][i]['nClusters'] for i in range(1,len(species_list))]

    # #species_diversity_cluster_distribution

    # import matplotlib.pyplot as plt
    # import numpy as np

    # hist, bins = np.histogram(species_diversity_cluster_distribution, bins=50)
    # width      = 0.7 * (bins[1] - bins[0])
    # center     = (bins[:-1] + bins[1:]) / 2

    # plt.bar(center, hist, align='center', width=width)
    # plt.ylabel('Frequency')
    # plt.xlabel('Universality (number of unique species)')
    # plt.title('Distribution of clusters by universality (e='+str(e)+', d='+str(save_d)+')')
    # plt.show()
