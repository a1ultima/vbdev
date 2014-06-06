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