# WARNINGS:
# CTRL+F "@WARN"

def load_motif_cluster_stats_dict():
    """
    Description: 

    Loads motif cluster statistics dictionary generated via: 
    scripts/stamp_to_collapsed/motifStatistics.p. See: http://goo.gl/1SQMyM

    Arguments:

    ...

    """
    import pickle
    f = open( '../../data/stamp_data/out/dreme_100bp_e0.05/SWU_SSD/cluster_to_stats.p' )
    cluster_to_stats = pickle.load( f )
    f.close()
    return cluster_to_stats

def get_motif_tree_collapsing_distances(cluster_to_stats):
    """
    Description:

    motif tree collapsed at a list of distances thresholds, see: scripts/stamp_to_collapsed/motifStatistics.p 

    Arguments:

    cluster_to_stats    see: load_motif_cluster_stats_dict()

    """
    parameterisation_events = cluster_to_stats.keys()
    distances = [float(i.split('d')[1]) for i in parameterisation_events]
    return distances

def blacklist_criteria_check(cluster, blacklist_motifs, cluster_to_stats, e, d, entropy_threshold=10.0, species_threshold=3 ):
    """
    Description:

    A motif_cluster is tested for CTRL+F "blacklisting criteria": 

    (1) Motifs of cluster are not present in the list of motifs of blacklisted clusters

    (2) cluster_H <= threshold_H     # low cluster entropy means "precise biological function"

    (3) cluster_S >= threshold_S     # high cluster species number (conservedness) means "less likely to be an artifact"

    Arguments:

    cluster = c412_n005_CCGAAYGCC   # string
    cluster_to_stats = load_motif_cluster_stats_dict() # 
    e = 0.05                        # string or float
    d = 0.16200000000000001         # string or float

    """

    # print cluster
    # cluster = 'c412_n005_CCGAAYGCC'
    # print type(cluster_to_stats)
    # print type(cluster_to_stats['e'+str(e)+'_d'+str(d)])
    # print type(cluster_to_stats['e'+str(e)+'_d'+str(d)][cluster])
    # print type(cluster_to_stats['e'+str(e)+'_d'+str(d)][cluster])

    # print type(cluster_to_stats['e'+str(e)+'_d'+str(d)][cluster]['motif'])
    # print type(cluster_to_stats['e'+str(e)+'_d'+str(d)][cluster]['motif']['list'])

    motifs = cluster_to_stats['e'+str(e)+'_d'+str(d)][cluster]['motif']['list'] # CTRL+F (1)


    species_number  = cluster_to_stats['e'+str(e)+'_d'+str(d)][cluster]['species']['unique']['n_unique'] # CTRL+F (2)

    entropy         = cluster_to_stats['e'+str(e)+'_d'+str(d)][cluster]['cluster']['H'] # CTRL+F (3)
 
    if not bool(set(motifs) & set(blacklist_motifs)):   # tests if any elements between "motifs" and "blacklist_motifs" overlap
                                                     # WARN: should test this, got it from StackOverflow: http://stackoverflow.com/a/3170067/3011648
                                                     # @TEST: try returning a separate label or print a warning whenever this condition is satisfied/not instead of just "False", because then we can discriminate it vs. the "False" that is supplied by the inner (2)/(3) criteria, e.g. maybe return the warning "blacklisted already" then diagnose the clusters
        print '\t\t\tCluster: '+cluster+' ...not blacklisted!'
        if ((entropy<=entropy_threshold) and (species_number>=species_threshold)):
            return True
        else:
            return False
    else:
        return False



#MAIN

print 'Loading pickled motif_cluster_statistics data...'
cluster_to_stats    = load_motif_cluster_stats_dict()
distances           = sorted(get_motif_tree_collapsing_distances(cluster_to_stats),reverse=True)


print 'Progressively gathering motif clusters satisfying entropy and species_number thresholds...'
# initiate...
blacklist           = [] # motifs satisfying the CTRL+F "blacklisting criteria"
blacklist_motifs    = [] # flat list of motifs of clusters blacklisted
                        # used to check if clusters of downstream distances
                        # are included in the upstream (lower distance,d) 
                        # clusters
                        # @WARN: this assumes all motifs of downstream d's cluster are included in upstream d's cluster?
                        # @TEST: 
                        #   <think of a test...> 

e = 0.05 # @E-VALUE parameter used by dreme to generate the motif data

entropy_threshold = 10.0

species_threshold = 3

for d in distances:

    print '\tdistance: '+str(d)

    clusters = cluster_to_stats['e'+str(e)+'_d'+str(d)].keys() # list of names of clusters, e.g. 

    for c in clusters:

        #print "\t\tcluster: "+str(c)

        if blacklist_criteria_check(c, blacklist_motifs, cluster_to_stats, e, d ,entropy_threshold,species_threshold):

            blacklist.append(c) # blacklist the cluster satisfying the 
                                # CTRL+F "blacklist critera"

            # Update the blacklist_motifs list, @WARN: test if the 
            motifs = cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['motif']['list']  # list of motifs belonging to cluster c
                    # e.g. cluster_to_stats['e0.05_d0.02']['c520_n027_GCTRTCA']['motif']['list']
                    # @WARN: all of this assumes an e0.05, but can easily be refactored to handle future e values, CTRL+F "@e-value"
            for m in motifs: 
                blacklist_motifs.append(m)

        else:
            pass

print blacklist










