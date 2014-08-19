# WARNINGS:
# CTRL+F "@WARN"
# TODOS:
# CTRL+F "@TODO"

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

    # @TODO: remove this debugging print stuff later
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
        # @TODO print junk
        #print '\t\t\tCluster: '+cluster+' ...not blacklisted!'
        if ((entropy<=entropy_threshold) and (species_number>=species_threshold)):
            print '\t\t\tCluster: '+cluster+' ...blacklisted!'
            return True
        else:
            return False
    else:
        return False


def get_blacklisted_motif_clusters( e = 0.05, entropy_threshold = 10.0, species_threshold = 3.0, cluster_to_stats=None):
    """

    Description:

    Returns a list, "blacklist", of motif clusters* that are canditates for trust-worthy de novo discovered motifs based on two criteria being satisfied: entropy lower than entropy_threshold, and species number lower than species_threshold. These criterial reflect the motif's biological functional precision and level of conservtion respectively. Conservation is used as a benchmark for filtering out motifs that are possibly the result of experimental errors. 

    *motif clusters are defined in /home/ab108/0VB/2kb/scripts/collapseMotifTree.py in exploreCutoffs()

    Arguments:

    e                   #   @E-VALUE parameter used by dreme to generate the motif data in <@TODO>

    entropy_threshold   #   see arguments for blacklist_criteria_check() <@TODO figure out a suitable value that will optimise for>

    species_threshold   #   see arguments for blacklist_criteria_check() <@TODO see above @TODO>

    cluster_to_stats    #   can be ignored, it is used as a debugging tool

    Returns:

    blacklist = [('c104_n004_CTGCGATCK', 0.5),
                 ('c106_n006_AGGGTAT', 0.5),
                 ('c029_n005_ACTGATGCA', 0.5),...] 

                 where 'c104_n004_CTGCGATCK' is the name of a motif cluster*, and 0.5 is the distance at which the motif tree* was collapsed (*see: /home/ab108/0VB/2kb/scripts/collapseMotifTree.py)

    cluster_to_stats    #   can be ignored, it is used as a debugging tool

    """

    print 'Loading pickled motif_cluster_statistics data...'
    if not cluster_to_stats:
        cluster_to_stats = load_motif_cluster_stats_dict()

    distances = sorted(get_motif_tree_collapsing_distances(cluster_to_stats),reverse=True)

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

    for d in distances:

        print '\tdistance: '+str(d)

        clusters = cluster_to_stats['e'+str(e)+'_d'+str(d)].keys() # list of names of clusters, e.g. 

        for c in clusters:

            #print "\t\tcluster: "+str(c)

            if blacklist_criteria_check(c, blacklist_motifs, cluster_to_stats, e, d ,entropy_threshold,species_threshold):

                # ENTANGLE c to d: cluster_to_stats unfortunately groups each cluster by parameters d and e. This grouping information is not implicit to the names of each cluster, e.g. c930_n005_AGTTGACGA. This means we cannot simply use the names of each cluster as keys to cluster_to_stats to retrieve informtion of that cluster; we need to associate the d of each cluster to the c of each cluster in order to retrieve information of that cluster. We can do this by placing c and it's d in a tuple, (c,d). These will be appended to the blacklist. 
                blacklist.append((c,d)) # blacklist the cluster satisfying the 
                                    # CTRL+F "blacklist critera"
                                    # Update the blacklist_motifs list
                motifs = cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['motif']['list']  # list of motifs belonging to cluster c
                        # e.g. cluster_to_stats['e0.05_d0.02']['c520_n027_GCTRTCA']['motif']['list']
                        # @WARN: all of this assumes an e0.05, but can easily be refactored to handle future e values, CTRL+F "@e-value"
                for m in motifs: 
                    blacklist_motifs.append(m)

            else:
                pass

    return blacklist, cluster_to_stats

def print_blacklist_summary(blacklist, cluster_to_stats, e=0.05):

    """

    Description:

    returns H_mean, S_mean: the mean entropy and mean species number of the motif clusters in the "blacklist", see get_blacklisted_motif_clusters().

    Arguments:

    blacklist = get_blacklisted_motif_clusters()

    cluster_to_stats = get_motif_tree_collapsing_distances() 

    e=0.05 # CTRL+F '@E-VALUE'

    """

    import numpy as np

    # MEAN SPECIES: Mean number of unique species in motif clusters of the blacklist
    S = [] # numbers of species per cluster in blacklist
    for (c,d) in blacklist: # (cluster,distance)
        s = cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['species']['unique']['n_unique'] # species number of cluster c
        S.append(s)
    S_mean = np.mean(S)


    # MEAN ENTROPY: Mean entropy of motif clusters of the blacklist
    H = []
    for (c,d) in blacklist:
        h = cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['cluster']['H'] # entropy of cluster c
        H.append(h)
    H_mean = np.mean(H)

    print '\n_________________________________'
    print 'SUMMARY STATISTICS OF BLACKLIST: '
    print '_________________________________\n'

    print '\tEntropy: '+str(H_mean)
    print '\tSpecies: '+str(S_mean)
    print '\tCluster: '+str(len(blacklist))+'         <-- no. of putative motifs'

    return H_mean, S_mean



e = 0.05

cluster_to_stats



blacklist, cluster_to_stats = get_blacklisted_motif_clusters( e = 0.05, entropy_threshold = 1.0, species_threshold = 3.0, cluster_to_stats=None) # --> 14 motifs, as expected 

summary = print_blacklist_summary(blacklist, cluster_to_stats, e=0.05


###

# @TESTS (get_blacklisted_motif_clusters):
#  
# if entropy_threshold is lower we should get <more/less> motifs?
#
# if species_threshold is higher we should get less motifs
#
#
# we get 325 motifs 

## ...
# blacklist, cluster_to_stats = get_blacklisted_motif_clusters( e = 0.05, entropy_threshold = 10.0, species_threshold = 3.0, cluster_to_stats=cluster_to_stats) # --> 325 motifs

## LOW SPECIES => LOW No. MOTIFS
# blacklist, cluster_to_stats = get_blacklisted_motif_clusters( e = 0.05, entropy_threshold = 10.0, species_threshold = 20.0, cluster_to_stats=cluster_to_stats) # --> 5 motifs, as expected

## LOW ENTROPY => LOW No. MOTIFS
#blacklist, cluster_to_stats = get_blacklisted_motif_clusters( e = 0.05, entropy_threshold = 1.0, species_threshold = 3.0, cluster_to_stats=cluster_to_stats) # --> 14 motifs, as expected

#blacklist, cluster_to_stats = get_blacklisted_motif_clusters( e = 0.05, entropy_threshold = 1.0, species_threshold = 3.0, cluster_to_stats=None) # --> 14 motifs, as expected

#summary = print_blacklist_summary(blacklist, cluster_to_stats, e=0.05)




# @TEST: try all combinations of input args to blacklist and test if it crashes it

###

# Summary Statistics for Bob:

# of a blacklist:
#   number of motifs
#   mean entropy
#   mean species number

import sys,getopt
    
# Argument Handling
def main(argv):

    # USER MADE ERROR? 

    # ... if no
    try:
        opts, args = getopt.getopt(argv,"hi:p:k:",["infile=","pvalue=","kcladetuple="])
    # ... if yes (e.g. mispelled args, non-existant args etc)
    except getopt.GetoptError:
        print 'fimo_report_expr_cladecombos.py -i <infile> -p <p-value threshold> -k <k number of clades>'
        sys.exit(2)

    # PROCESS THE @ARGS

    for opt, arg in opts:
        if opt == '-h': # help file
            print 'fimo_report_expr_cladecombos.py -i <infile> -p <p-value threshold> -k <k number of clades>'
            sys.exit()
        elif opt in ("-i","--infile"):    # long and short args, see
            i = arg
        elif opt in ("-p", "--pvalue"):
            p = arg
        elif opt in ("-k", "--kcladetuple"):
            k = arg

    # CMD or PYTHON EXECUTE CONTROL: by default try run as if on cmd line, with input and output args ... else if it errors it may be running from within python, so then run default args: 

    # ... script was run from command line
    try:
        
        # SCRIPT
         # <-- main script, @REPLICATE exactly as below
        #                             ^@args here processed from above
    
    # ... script was run within Python
    except UnboundLocalError:
        
        i = raw_input('input file? e.g. /home/maccallr/agcc/fimo/merged-motifs-d0.002-species3/all-concatenated/fimo_anopheles_gambiae/fimo.report-expr.txt')
        
        p = raw_input('p-value threshold? e.g. 0.1')
        
        k = raw_input('How many k clade tuples (0-3)? e.g. 2')

        fimo_report_expr_cladecombos(i,p,k)  # <-- main script, CTRL+F "@replicate" exactly as above

if __name__ == "__main__":
    main(sys.argv[1:])




