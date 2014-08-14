
# MOTIF STATS
import os 
print(os.getcwd())
if not os.getcwd().endswith('scripts'):
    os.chdir('../')
    print os.getcwd()

import imp
speciesManage = imp.load_source('speciesManage.py', './speciesManage.py')
motif_entropy = imp.load_source('motif_entropy.py', './motif_entropy.py')

import numpy as np 
import copy

#import matplotlib.pyplot as plt

event_to_clusterClass   = {}                        # events as primary key


# Ask user to provide the e-value cut-off pointing to the data we need, e.g. 0.5 points to dreme_100bp_e0.5. Then we get a list of d-cutoff values available in the directory by re-formatting the filenames
e_cut       = str(raw_input('what is the e-value of the data you want to work with? e.g. 0.5'))
e_cut_dir   = '../data/stamp_data/out/dreme_100bp_e'+e_cut+'/SWU_SSD/'
d_list      = [float(i.replace('cluster_motifs_d','')) for i in os.listdir(e_cut_dir) if 'cluster_motifs_d' in i]


parameterisation_events = zip([float(e_cut)]*len(d_list),d_list) # list of 'parameterisation events', e.g. 'of a parameterisation event': ( evalue, dvalue ), real examples commented below
#parameterisation_events = [(0.05,0.05),(0.05,0.5)]     
#parameterisation_events = [(0.05,0.5)]                 

filename_species= './species_list.txt'
species_list    = speciesManage.generate_list(filename_species) # generates a list of species names corresponding to EnsEMBl MySQL db name

species_to_index= {} # Speices:Index dict: species now point to the indices of a list which will later be the vector of counts and presence/absence for phylo profiles

for i,s in enumerate(species_list):
    species_to_index[s] = i

# 
# MAIN LOOP             :   iteratre through event combinations in 'parameterisation_events', e.g. [] evalue
#

cluster_to_stats = {} 

for e,save_d in parameterisation_events:

    event = 'e'+str(e)+'_d'+str(save_d)

    print('Parameters: e:'+str(e)+' d:'+str(save_d))

    # symbol_to_species dict    : to later translate motifs names into species name (CTRL+F: SPECIES-to-MOTIFS  ) => to generate dict of species_to_motifs => to then 
                                # e.g. 
                                #   of a symbol: DMEL_, 
                                #   of a species name: Drosophila melanogaster
    symbol_to_species = {}
    for species in species_list:
        symbol = species.split('_')[0][0].upper()+species.split('_')[1][0:3].upper()
        symbol_to_species[symbol] = species

    #------------------------------------------------------------------------------
    # SPECIES-to-MOTIFS 
    #
    # species_to_motifs         : to get #motifs per species
    filename_dreme      = '../data/stamp_data/in/dreme_100bp_e'+str(e)+'/dreme.fasta'
    file_dreme          = open(filename_dreme,'r') 
    
    species_to_motifs   = {} # species as key to motif list

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

    print 'motifs'

    species_nmotifs_dist = [] # species vs. number of motifs distribution
    
    for species in species_to_motifs.keys(): # iterate through species
        nmotifs = len(species_to_motifs[species])   # count number of motifs on given species
        species_nmotifs_dist.append(nmotifs)        # append to the species vs. number of motifs distribution
        print('\t'+species + ': ' + str(len(species_to_motifs[species])))
    
    species_nmotifs_dist = np.array(nmotifs)

#------------------------------------------------------------------------------------
#
# MOTIF CLUSTER STATISTICS:
#
#   Motif cluster statistics grouped/contextualised in two main ways:
# 
#   (1) Cluster statistics represented in 'parameterisation event':'cluster':<x...> where the list of X ways of grouping (classes) is found by CTRL+F: "cluster_to_stats[event][cluster]"
#
#   (2) Cluster statistics represented in n-species-centric manner, where n represents a class: the number of species OR MORE that motif clusters are grouped by, e.g. "motifs clusters of which there are 3 or more species"
#

    #------------------------------------------------------------------------------
    #  (1) ...
    #
    #   Generate various statistics for the clusters, in various contexts, e.g. in the context of their species, e.g. species 1: <motif cluster statistics>, species 2: <motif cluster statistics>, ...
    #
    #   CTRL+F: "cluster_to_stats[event][cluster]" for more specific idea of how the data is structured
    # 
    clusters_fbp_path       = '../data/stamp_data/out/dreme_100bp_e'+str(e)+'/SWU_SSD/cluster_motifs_d'+str(save_d)+'/family_binding_matrices/'
    clusters_tree_path      = '../data/stamp_data/out/dreme_100bp_e'+str(e)+'/SWU_SSD/cluster_motifs_d'+str(save_d)+'/family_binding_matrices/'
    clusters_jaspar_path    = '../data/stamp_data/out/dreme_100bp_e'+str(e)+'/SWU_SSD/cluster_motifs_d'+str(save_d)+'/family_binding_matrices/'
    clusters_transfac_path  = '../data/stamp_data/out/dreme_100bp_e'+str(e)+'/SWU_SSD/cluster_motifs_d'+str(save_d)+'/family_binding_matrices/'
    clusters_fasta_path     = '../data/stamp_data/out/dreme_100bp_e'+str(e)+'/SWU_SSD/cluster_motifs_d'+str(save_d)+'/fasta/'
    clusters_path           = '../data/stamp_data/out/dreme_100bp_e'+str(e)+'/SWU_SSD/cluster_motifs_d'+str(save_d)+'/headers_for_fasta/'

    clusters_files          = os.listdir(clusters_path)
    species_to_clusters     = {}
    cluster_to_stats[event] = {}
    cluster_entropies       = []

    for count,cluster_filename in enumerate(clusters_files):
        
        cluster = cluster_filename.replace('.txt','')    
        
        # 
        # 
        #
        cluster_to_stats[event][cluster] = { 
            'species':{
                'unique':{
                    'n_unique':0,
                    'n_paralogue':0,
                    'list':[]
                    },
                'phy_profile':{
                    'binary':{'dict_mode':{},'list_mode':[0]*len(species_list)},
                    'count':{'dict_mode':{},'list_mode':[0]*len(species_list)}
                    }
                },
            'motif':{
                'n':0,
                'list':[],
                'avg_length':[]
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
        for count_sp,species in enumerate(species_list):
            cluster_to_stats[event][cluster]['species']['phy_profile']['binary']['dict_mode'][species] = 0
            cluster_to_stats[event][cluster]['species']['phy_profile']['count']['dict_mode'][species]  = 0

    # Initiate filenames
        # Filenames  
        cluster_tree    = clusters_tree_path+cluster+'.tree'
        cluster_jaspar  = clusters_jaspar_path+cluster+'_match_paris.txt'
        cluster_transfac= clusters_transfac_path+cluster+'_matched.transfac'
    
        # If the cluster has only one motif, then it does not have an FBP.txt file, so seek for it's original TRANSFAC fasta file to represent the FBP
        if 'n001' in cluster:   # ... FBP.txt
            cluster_fbp                                 = clusters_fasta_path+cluster+'.fasta'
            cluster_fbp_entropy                         = motif_entropy.fbp(cluster_fbp,format='count')
            cluster_to_stats[event][cluster]['cluster']['H']   = cluster_fbp_entropy # calculate entropy for the motif
        else:                   # ... .fasta
            cluster_fbp                                 = clusters_fbp_path+cluster+'FBP.txt'
            cluster_fbp_entropy                         = motif_entropy.fbp(cluster_fbp,format='decimal')
            cluster_to_stats[event][cluster]['cluster']['H']   = cluster_fbp_entropy # calculate entropy for the motif

        cluster_to_stats[event][cluster]['cluster']['path']['fbp']     = cluster_fbp
        cluster_to_stats[event][cluster]['cluster']['path']['tree']    = cluster_tree
        cluster_to_stats[event][cluster]['cluster']['path']['jaspar']  = cluster_jaspar
        cluster_to_stats[event][cluster]['cluster']['path']['transfac']= cluster_transfac

    # CLUSTER HEADERS           : crawl through cluster file to read the motif names
        cluster_file    = open(clusters_path+cluster_filename,'r')
        cluster_species = []
        motif_lenghts   = []

        while True: # iterates through each motif in the cluster
            motif = cluster_file.readline().rstrip()
            if motif == "":
                break
            symbol  = motif.split('_')[0]
            try: 
                species = symbol_to_species[symbol]
            except KeyError:
                raise Exception('There is a missing Key, check that the species_list.txt matches with the species available in the data for this event... ')
            
            cluster_species.append(species)
            
            cluster_to_stats[event][cluster]['motif']['n'] = 1 + cluster_to_stats[event][cluster]['motif']['n']
            
            cluster_to_stats[event][cluster]['motif']['list'].append(motif)
            
            motif_lenghts.append(len(motif.split('_')[1])) # collect length of the current motif to later average them per cluster
        cluster_file.close()

    # AVERAGE MOTIF LENGHT
        cluster_to_stats[event][cluster]['motif']['avg_length'] = np.array(motif_lenghts).mean() # average the lengths of motifs

    # append the cluster to species who were found in the cluster
        for count_sp,species in enumerate(cluster_species):
            cluster_to_stats[event][cluster]['species']['phy_profile']['binary']['dict_mode'][species] = 1
            cluster_to_stats[event][cluster]['species']['phy_profile']['count']['dict_mode'][species]  = 1 + cluster_to_stats[event][cluster]['species']['phy_profile']['count']['dict_mode'][species]

        for s in cluster_to_stats[event][cluster]['species']['phy_profile']['binary']['dict_mode'].keys():
            index   = species_to_index[s]
            binary  = cluster_to_stats[event][cluster]['species']['phy_profile']['binary']['dict_mode'][s] # presence/absence of species s in the dict
            count1  = cluster_to_stats[event][cluster]['species']['phy_profile']['count']['dict_mode'][s] # counts of ^^ of species 

            cluster_to_stats[event][cluster]['species']['phy_profile']['binary']['list_mode'][index]   = binary
            cluster_to_stats[event][cluster]['species']['phy_profile']['count']['list_mode'][index]    = count1

        unique_cluster_species = set(cluster_species)    # get the unique list of species within cluster

        # cluster_to_species
        cluster_to_stats[event][cluster]['species']['unique']['n_unique']      = len(unique_cluster_species)
        cluster_to_stats[event][cluster]['species']['unique']['n_paralogue']   = len(cluster_species) - len(unique_cluster_species)
        cluster_to_stats[event][cluster]['species']['unique']['list']          = unique_cluster_species

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

        # @PRINT
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

        event_to_clusterClass[event] = copy.deepcopy(event_to_clusterClass_template) # cluster classes

        # all               :   summary stats for all the clusters
        clusters_all                                                = cluster_to_stats[event].keys()
        
        clusters_all_entropy                                        = np.array([cluster_to_stats[event][c]['cluster']['H'] for c in clusters_all]) # entropies for each cluster
        event_to_clusterClass[event]['all']                         = copy.deepcopy(clusterClass_to_stats_template)
        event_to_clusterClass[event]['all']['list']                 = clusters_all
        event_to_clusterClass[event]['all']['nClusters']            = len(clusters_all)
        event_to_clusterClass[event]['all']['H']['distribution']    = clusters_all_entropy
        event_to_clusterClass[event]['all']['H']['mean']            = clusters_all_entropy.mean()

        # nSpeciesOrMore    :   clusters w/ species >= X
        for s in range(1,len(species_list)+1): # per species number, S
            
            clusters_s                                                              = [c1 for c1 in cluster_to_stats[event].keys() if cluster_to_stats[event][c1]['species']['unique']['n_unique'] >= s] # find clusters with nSpecies >= s
            clusters_s_entropy                                                      = np.array([cluster_to_stats[event][c]['cluster']['H'] for c in clusters_s]) # entropies for each cluster

            event_to_clusterClass[event]['nSpeciesOrMore'][s]                       = copy.deepcopy(clusterClass_to_stats_template)
            event_to_clusterClass[event]['nSpeciesOrMore'][s]['list']               = clusters_s
            event_to_clusterClass[event]['nSpeciesOrMore'][s]['nClusters']          = len(clusters_s)
            event_to_clusterClass[event]['nSpeciesOrMore'][s]['H']['distribution']  = clusters_s_entropy
            event_to_clusterClass[event]['nSpeciesOrMore'][s]['H']['mean']          = clusters_s_entropy.mean()

        print(len(event_to_clusterClass[event]['all']['list']))

        # clades            :   one level deeper than the above, under clades we have various clade related conditions to define the classes the template must be appended to each
        for clade in event_to_clusterClass_template['clades'].keys():
            event_to_clusterClass[event]['clades'][clade] = copy.deepcopy(clusterClass_to_stats_template) # append common data template to each clade

        # inc_DMEL          :   clusters including at least one drosophila
        clusters_inc_DMEL                                                       = [c for c in cluster_to_stats[event].keys() if 'droso' in ''.join(cluster_to_stats[event][c]['species']['unique']['list'])] # find clusters that include drosophila motif/s

        clusters_inc_DMEL_entropy                                               = np.array([cluster_to_stats[event][c]['cluster']['H'] for c in clusters_inc_DMEL]) # entropies for each cluster
        event_to_clusterClass[event]['clades']['inc_DMEL']['list']              = clusters_inc_DMEL
        event_to_clusterClass[event]['clades']['inc_DMEL']['nClusters']         = len(clusters_inc_DMEL)
        event_to_clusterClass[event]['clades']['inc_DMEL']['H']['distribution'] = clusters_inc_DMEL_entropy
        event_to_clusterClass[event]['clades']['inc_DMEL']['H']['mean']         = clusters_inc_DMEL_entropy.mean()

        #clusters_inc_gambiaeComplex = 

    #------------------------------#
    #   Data for the SpreadSheet   #
    #------------------------------#

    # CLUSTERS:

    # Names of clusters (column 1)
    cluster_names       = cluster_to_stats[event].keys()

    # Numbers of motifs (column 2)
    cluster_nMotifs     = [cluster_to_stats[event][i]['motif']['n'] for i in cluster_to_stats[event].keys()]

    # Entropy of clusters ()
    cluster_entropies   = [cluster_to_stats[event][i]['cluster']['H'] for i in cluster_to_stats[event].keys()]

    # Numbers of unique species ()
    cluster_nSpecies   = [cluster_to_stats[event][i]['species']['unique']['n_unique'] for i in cluster_to_stats[event].keys()]

    # Numbers of paralogues ()
    cluster_nParalogues= [cluster_to_stats[event][i]['species']['unique']['n_paralogue'] for i in cluster_to_stats[event].keys()]

    #   PHYLOGENETIC PROFILE:

    # Phylo profile - binary ()
    cluster_phy_binary= [cluster_to_stats[event][i]['species']['phy_profile']['binary']['list_mode'] for i in cluster_to_stats[event].keys()]

    # Phylo profile - count ()
    cluster_phy_count= [cluster_to_stats[event][i]['species']['phy_profile']['count']['list_mode'] for i in cluster_to_stats[event].keys()]

    # EVENTS:
    #   CLUSTERS >=S: 

    # Entropies for clusters >= S ()
    cluster_sMore_meanEntropies = [event_to_clusterClass[event]['nSpeciesOrMore'][s]['H']['mean'] for s in event_to_clusterClass[event]['nSpeciesOrMore'].keys()] # when you plot entropy vs. nSpecies the points higher than the diagonal are those which are high confidence motifs, i.e. low entropy but present in many species 

    # Number of motifs for clusters >= S ()
    cluster_sMore_nMotifs = [event_to_clusterClass[event]['nSpeciesOrMore'][s]['nClusters'] for s in event_to_clusterClass[event]['nSpeciesOrMore'].keys()] 

    #   HOLLISTIC:

    # Total number of clusters
    all_nClusters = event_to_clusterClass[event]['all']['nClusters']

    # Average entropy of clusters
    all_avgEntropy = event_to_clusterClass[event]['all']['H']['mean']

    # Average entropy of clusters
    all_avgEntropy = event_to_clusterClass[event]['all']['H']['mean']


# Save cluster_to_stats[event][cluster] dict to: /home/ab108/0VB/2kb/data/stamp_data/out/dreme_100bp_e0.05/SWU_SSD/cluster_to_stats.p
import pickle
pickle_dir = e_cut_dir+'/cluster_to_stats'
pickle_obj = cluster_to_stats
pickle.dump(pickle_obj,open(pickle_dir+'.p','wb'))
