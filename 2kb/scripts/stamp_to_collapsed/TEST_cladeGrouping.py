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

    blacklist = get_blacklisted_motif_clusters()                # 

    cluster_to_stats = get_motif_tree_collapsing_distances()    #

    e=0.05 # CTRL+F '@E-VALUE' 

    """

    import numpy as np

    # MEAN SPECIES: Mean number of unique species in motif clusters of the blacklist
    S = [] # numbers of species per cluster in blacklist

    for (c,d) in blacklist: # (/cluster,distance)

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
    print '\tCluster: '+str(len(blacklist))+'         <-- no. of putative motifs\n\n\n'

    return H_mean, S_mean



def blacklist_then_summaryStats_from_shell( e = 0.05, species_threshold = 3.0, entropy_threshold = 1.0 ):
    """ 

    Description:\n\n

    A wrapper for obtaining the blacklist at given species number threshold, S, and entropy threshold, H, (CTRL+F: blacklist_criteria_check) and e value (CTRL+F: @e-value)\n\n

    Arguments:\n\n 

    e = 0.05                      # CTRL+F: @e-value\n\n

    species_threshold = 3.0       # CTRL+F: blacklist_criteria_check\n\n

    entropy_threshold = 1.0       # CTRL+F: blacklist_criteria_check\n\n

    Example:\n\n

    blacklist_then_summaryStats_from_shell( e = 0.05, species_threshold=3.0, entropy_threshold=1.0)\n\n

    """

    blacklist, cluster_to_stats = get_blacklisted_motif_clusters( e, species_threshold, entropy_threshold, cluster_to_stats=None) # --> 14 motifs, as expected 

    summary = print_blacklist_summary( blacklist, cluster_to_stats, e)

    return blacklist, summary



#----------------------------------------

#------------------------------------------------------------------------











"""

blacklisted_cluster = blacklist[6] # ('c087_n073_ACGCACGMA', 0.5) is the most species rich cluster


max_len = 0

for j,i in enumerate(blacklist):

    c=i[0]
    d=i[1]
    unique_species = cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['species']['unique']['list'] # e.g. ['anopheles_maculatus','anopheles_atroparvus','anopheles_sinensis','anopheles_melas','anopheles_merus','anopheles_epiroticus','aedes_aegypti','anopheles_darlingi','anopheles_gambiae','anopheles_culicifacies','anopheles_dirus','culex_quinquefasciatus','anopheles_arabiensis','anopheles_farauti','anopheles_quadriannulatus','anopheles_funestus','anopheles_minimus','anopheles_christyi','anopheles_stephensi','anopheles_albimanus']

    tmp_len = len(unique_species)

    if tmp_len > max_len:
        max_len = tmp_len
        species_rich_cluster = (c,d)
        index_sp = j

    # if 'drosophila_melanogaster' in unique_species: # @TESTED: there are drosophila clusters!
    #     print i
    # else:
    #     print 'lol'
"""


#---------


# MAKE THE IF STATEMENTS FOR CLADES

# species from cluster to stats 


def blacklist_to_clades( blacklist ):

    """

    Description: 

    Returns five lists of @motif clusters* grouped by clade**. Each list is a group of motif clusters based on the branch of the phylogenetic tree that their motifs were originally predicted in the @promoterome (s) by dreme.     

    * motif clusters are in the format (clustername, distance***), e.g. ('c087_n073_ACGCACGMA', 0.5)

    ** clade-based groups of clusters are: "dipteran_clusters", "mosquito_clusters", "anopheles_clusters" and "gambiae_clusters":

            dipteran_clusters: contains clusters if cluster has at least one motif belonging to a dipteran (drosophila), i.e. very conserved motifs, i.e.e. evolved early on in dipteran ancestry. see: cluster_has_dipteran().

            mosquito_clusters: clusters having at least one motif belonging to a mosquito species (aedes, culex, others, but not other dipterans). These are less conserved than dipteran clusters. see: cluster_has_mosquito().

            anopheles_clusters: clusters having at least one motif belonging to an anopheline species plus others, but not other dipterans or mosquitoes (for the list see: cluster_has_anophelines() ). see: cluster_has_anophelines().

            gambiae_clusters: clusters having at least one motif belonging to gambiae complex species, but not other dipterans or mosquitoes or anophelines. see: cluster_has_gambiae_complex().

            Note: these groups are independent, i.e. no overlaps in clusters. See descriptions of each groups personalised sub-functions for more info. 


    *** see CTRL+F "0.5 is the distance at" to understand what we mean by distance.


    Arguments: 

    blacklist, ... = get_blacklisted_motif_clusters( e, entropy_threshold, species_threshold, cluster_to_stats=None)  # CTRL+F: "@blacklist"


    """

    #------------------------------------------------------------------------
    # CLADE CHECKER FUNCTIONS   

    def cluster_has_dipteran( blacklisted_cluster ):

        """

        Description: 

        returns TRUE if the blacklisted_cluster contains at least one Dipteran* @motif**

        * Dipteran:
            drosophila_melanogaster 

        ** @MOTIF: a motif that was originally predicted by dreme when fed the drosophila @promoterome** as input data

        *** @PROMOTEROME: is the sequence data generated by 2kb/scripts/upDownStream.py


        """

        c = blacklisted_cluster[0] # cluster name
        d = blacklisted_cluster[1] # distance where it came from

        unique_species = cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['species']['unique']['list'] # e.g. ['anopheles_maculatus','anopheles_atroparvus','anopheles_sinensis','anopheles_melas','anopheles_merus','anopheles_epiroticus','aedes_aegypti','anopheles_darlingi','anopheles_gambiae','anopheles_culicifacies','anopheles_dirus','culex_quinquefasciatus','anopheles_arabiensis','anopheles_farauti','anopheles_quadriannulatus','anopheles_funestus','anopheles_minimus','anopheles_christyi','anopheles_stephensi','anopheles_albimanus']

        if 'drosophila_melanogaster' in unique_species: # @TESTED: there are drosophila clusters, and this does find them
            return True
        else: 
            return False  

    def cluster_has_mosquito( blacklisted_cluster ):
        """

        Description: 

        returns TRUE if the blacklisted_cluster contains at least one Mosquito* @motif

        * Mosquito:
            culex_quinquefasciatus
            aedes_aegypti

        ** a motif that was originally predicted by dreme when fed the drosophila @PROMOTEROME** as input data

        *** @PROMOTEROME: is the sequence data generated by 2kb/scripts/upDownStream.py

        Arguments:

        @cluster
        
        """

        c = blacklisted_cluster[0] # cluster name
        d = blacklisted_cluster[1] # distance that motif-tree was collapsed to get the cluster

        unique_species = cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['species']['unique']['list'] # e.g. ['anopheles_maculatus','anopheles_atroparvus','anopheles_sinensis','anopheles_melas','anopheles_merus','anopheles_epiroticus','aedes_aegypti','anopheles_darlingi','anopheles_gambiae','anopheles_culicifacies','anopheles_dirus','culex_quinquefasciatus','anopheles_arabiensis','anopheles_farauti','anopheles_quadriannulatus','anopheles_funestus','anopheles_minimus','anopheles_christyi','anopheles_stephensi','anopheles_albimanus']

        if ('culex_quinquefasciatus' in unique_species) or ('aedes_aegypti' in unique_species): # @TESTED: there are drosophila clusters, and this does find them
            return True
            #print 'true'
        else: 
            return False
            #print 'false'

        # @TEST                     PASSED
        # if cluster has:
        #   only culex
        #       unique_species = ['culex_quinquefasciatus']
        #   only aedes
        #       unique_species = ['aedes_aegypti']
        #   neither
        #       unique_species = ['drosophila']
        #   both
        #       unique_species = ['culex_quinquefasciatus','aedes_aegypti']


    def cluster_has_anophelines( blacklisted_cluster ):

        """

        Description: 

        returns TRUE if the blacklisted_cluster contains at least one Anopheline* @motif

        Anopheline:
            anopheles_darlingi
            anopheles_albimanus
            anopheles_atroparvus
            anopheles_sinensis
            anopheles_farauti
            anopheles_dirus
            anopheles_minimus
            anopheles_funestus
            anopheles_culicifacies
            anopheles_maculatus
            anopheles_stephensi
            anopheles_epiroticus
            anopheles_christyi

        Arguments:

        @cluster

        """

        c = blacklisted_cluster[0] # cluster name
        d = blacklisted_cluster[1] # distance that motif-tree was collapsed to get the cluster

        unique_species = cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['species']['unique']['list'] # e.g. ['anopheles_maculatus','anopheles_atroparvus','anopheles_sinensis','anopheles_melas','anopheles_merus','anopheles_epiroticus','aedes_aegypti','anopheles_darlingi','anopheles_gambiae','anopheles_culicifacies','anopheles_dirus','culex_quinquefasciatus','anopheles_arabiensis','anopheles_farauti','anopheles_quadriannulatus','anopheles_funestus','anopheles_minimus','anopheles_christyi','anopheles_stephensi','anopheles_albimanus']

        anophelines_not_including_gambiae_complex = ['anopheles_darlingi','anopheles_albimanus','anopheles_atroparvus','anopheles_sinensis','anopheles_farauti','anopheles_dirus','anopheles_minimus','anopheles_funestus','anopheles_culicifacies','anopheles_maculatus','anopheles_stephensi','anopheles_epiroticus','anopheles_christyi'] # @TEST: check this list is correct, i.e. (i)contains all the anophelines, (ii) does not have any non-anophelines, (iii) does not have any gambiae complex

        is_overlaps_between_lists = bool(set(anophelines_not_including_gambiae_complex) & set(unique_species)) # checks for overlaps between the two lists, i.e. if the clusters species list contains any anophelines

        if is_overlaps_between_lists: # @TESTED: there are drosophila clusters, and this does find them
            return True
        else: 
            return False


    def cluster_has_gambiae_complex( blacklisted_cluster ):

        """
        Description: 

        returns TRUE if the blacklisted_cluster contains at least one Gambiae complex* @motif

        * Gambiae complex:
            anopheles_melas
            anopheles_merus
            anopheles_quadriannulatus
            anopheles_arabiensis
            anopheles_gambiae

        Arguments:

        @cluster

        """

        c = blacklisted_cluster[0] # cluster name
        d = blacklisted_cluster[1] # distance that motif-tree was collapsed to get the cluster

        unique_species = cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['species']['unique']['list'] # e.g. ['anopheles_maculatus','anopheles_atroparvus','anopheles_sinensis','anopheles_melas','anopheles_merus','anopheles_epiroticus','aedes_aegypti','anopheles_darlingi','anopheles_gambiae','anopheles_culicifacies','anopheles_dirus','culex_quinquefasciatus','anopheles_arabiensis','anopheles_farauti','anopheles_quadriannulatus','anopheles_funestus','anopheles_minimus','anopheles_christyi','anopheles_stephensi','anopheles_albimanus']

        gambiae_complex_not_other_anophelines = ['anopheles_melas','anopheles_merus','anopheles_quadriannulatus','anopheles_arabiensis','anopheles_gambiae']

         # @TEST: check this list is correct, i.e. (i)contains all the gambiae complex, (ii) does not have any other anophelines, (iii) does not have any gambiae complex

        is_overlaps_between_lists = bool(set(gambiae_complex_not_other_anophelines) & set(unique_species)) # checks for overlaps between the two lists, i.e. if the clusters species list contains any anophelines

        if is_overlaps_between_lists: # @TESTED: there are drosophila clusters, and this does find them
            return True
        else: 
            return False




    #---------



    # INITIATE CLADE GROUPINGS

    dipteran_clusters = [] # contains at least one drosophila 
    mosquito_clusters = [] # (contains no drosophila) AND (contains at least one culex OR one aedes)
    anopheles_clusters = [] 
    # (contains no mosquito AND no dipteran) AND (at least one anopheles (?) ) @BOB: wouldn't this just swallow up the whole of the gambiae complex...?
    gambiae_clusters = [] # contains (no anopheles AND no dipteran AND no mosquito) AND (at least one gambiae complex)

    # PER BLACKLISTED CLUSTER...

    for c in blacklist:

        # GROUP THE CLUSTER INTO ONE OF THE CLADEs

        if cluster_has_dipteran(c): # store cluster if it contains at least one dipteran (drosophila)
            dipteran_clusters.append(c)

        if ( cluster_has_mosquito(c) and (not cluster_has_dipteran(c)) ): # store cluster if it contains at least one mosquito, but not if it has a dipteran

        # @TEST: any difference in "and not" vs. "and ( not ... "  ??
            mosquito_clusters.append(c)

        if cluster_has_anophelines(c) and ( (not cluster_has_dipteran(c)) and (not cluster_has_mosquito(c))):  # store cluster if it contains at least one anopheline, but not if it has a dipteran or a mosquito

            anopheles_clusters.append(c)

        if cluster_has_gambiae_complex(c) and ( (not cluster_has_dipteran(c)) and (not cluster_has_mosquito(c)) and (not cluster_has_anophelines(c))): # store cluster if it contains at least one gambiae complex, but not if it has an anopheline or mosquito or dipteran
        #@TESTs_1: 
        #  - if the remainder == the output of this test, if yes, then great!
        #  - are there any overlaps between clusters in each of the groups? Hint: use set() 
            gambiae_clusters.append(c)

    # dipteran_clusters # contains at least one drosophila 

    # mosquito_clusters # (contains no drosophila) AND (contains at least one culex OR one aedes)

    # anopheles_clusters # (contains no mosquito AND no dipteran) AND (at least one anopheles (?) ) @BOB: wouldn't this just swallow up the whole of the gambiae complex...?

    # gambiae_clusters

    print '\n_________________________________'
    print 'CLADE DISTRIBUTION OF BLACKLIST: '
    print '_________________________________\n'
    print '\tDipteran: '    +str(len(dipteran_clusters))
    print '\tMosquito: '    +str(len(mosquito_clusters))
    print '\tAnopheles: '   +str(len(anopheles_clusters))
    print '\tGambiae: '     +str(len(gambiae_clusters))

    return dipteran_clusters, mosquito_clusters, anopheles_clusters, gambiae_clusters


#----------------------------------------



# TESTING for blacklist_to_clades(), see: CTRL+F "@TEST_1"
#
#   - tests passed for: e = 0.05, species_threshold = 3.0, entropy_threshold = 10.0
#
"""
e = 0.05                      # CTRL+F: @e-value\n\n
species_threshold = 3.0       # CTRL+F: blacklist_criteria_check\n\n
entropy_threshold = 10.0       # CTRL+F: blacklist_criteria_check\n\n
blacklist, cluster_to_stats = get_blacklisted_motif_clusters( e, entropy_threshold, species_threshold, cluster_to_stats=None) # If errors change: cluster_to_stats=None 
dipteran_clusters, mosquito_clusters, anopheles_clusters, gambiae_clusters = blacklist_to_clades(blacklist)
import itertools
groups = (dipteran_clusters, mosquito_clusters, anopheles_clusters, gambiae_clusters)
pairs = list(itertools.product(groups,repeat=2))
for i,p in enumerate(pairs):
    if bool( set(p[0]) & set(p[1]) ):   # if there are overlaps...
        if not p[0]==p[1]:              # ...AND the overlapping lists compared are not (self,self) [i.e. the diagonal of the 2-tuples matrix]
            print 'uh oh'               # ...then we are worried, it means our clade groups overlap in collapsed motif clusters, it means there is a silent bug in blacklist_to_clades()
"""


e = 0.05                      # CTRL+F: @e-value\n\n
species_threshold = 3.0       # CTRL+F: blacklist_criteria_check\n\n
entropy_threshold = 10.0       # CTRL+F: blacklist_criteria_check\n\n

blacklist, cluster_to_stats = get_blacklisted_motif_clusters( e, entropy_threshold, species_threshold, cluster_to_stats=cluster_to_stats) # If errors change: cluster_to_stats=None 


moo = blacklist_to_clades( blacklist )
