"""

Description:

...

Run the following scripts in order, in the way shown:

1. collapseMotifTree.py 
2. motifStatistics.py
3. python collapseMotifTree_progressiveMode.py --nspecies 3.0 --entropy 10.0

Notes:

Some technical words may not make sense to you. Such words may have intuitive descriptions available in here. Words with intuitive descriptions available will be in the format: @<word-in-lowercase>. Simply CTRL+F "@<word-in-uppercase>" to jump to the intuitive descriptions.


@standard examples:

cluster_to_stats['e0.05_d0.056']['c346_n002_CAATACWCG']

                    ^e = 0.05
                          ^d = 0.056  ^c = 'c346_n002_CAATACWCG'

"""

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

    print 'Loading motif cluster data structure, please wait a couple of minuites...'

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

    # ACCESS DATA

    motifs = cluster_to_stats['e'+str(e)+'_d'+str(d)][cluster]['motif']['list'] # CTRL+F (1)

    species_number  = cluster_to_stats['e'+str(e)+'_d'+str(d)][cluster]['species']['unique']['n_unique'] # CTRL+F (2)

    entropy         = cluster_to_stats['e'+str(e)+'_d'+str(d)][cluster]['cluster']['H'] # CTRL+F (3)
 

    # CHECK IF PARENT DISTANCE ALREADY BLACKLISTED

    if not bool(set(motifs) & set(blacklist_motifs)):   # tests if any elements between "motifs" and "blacklist_motifs" overlap
                                                     # @WARN: should test this, got it from StackOverflow: http://stackoverflow.com/a/3170067/3011648
                                                     # @TEST: try returning a separate label or print a warning whenever this condition is satisfied/not instead of just "False", because then we can discriminate it vs. the "False" that is supplied by the inner (2)/(3) criteria, e.g. maybe return the warning "blacklisted already" then diagnose the clusters

        # IS THE CLUSTER WORTHY? ...
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

    Returns two pieces of data: 
        (i) "@blacklist"* of motif @clusters** that are canditates for trust-worthy de novo discovered motifs based on two criteria being satisfied: entropy lower than entropy_threshold, and species number lower than species_threshold. These criterial reflect the motif's biological functional precision and level of conservtion respectively. Conservation is used as a benchmark for filtering out motifs that are possibly the result of experimental errors. 
        (ii) cluster_to_stats (ignore, it's a debugging tool)

    * blacklist = [  ('c104_n004_CTGCGATCK', 0.5),
                     ('c106_n006_AGGGTAT', 0.5),
                     ('c029_n005_ACTGATGCA', 0.5),...] 

                 where 'c104_n004_CTGCGATCK' is the name of a motif cluster*, and 0.5 is the distance at which the motif tree* was collapsed (*see: /home/ab108/0VB/2kb/scripts/collapseMotifTree.py)

    ** @CLUSTER: are defined in /home/ab108/0VB/2kb/scripts/collapseMotifTree.py in exploreCutoffs()
  

    Arguments:

    e                   #   @E-VALUE parameter used by dreme to generate the motif data in <@TODO>

    entropy_threshold   #   see arguments for blacklist_criteria_check() <@TODO figure out a suitable value that will optimise for>

    species_threshold   #   see arguments for blacklist_criteria_check() <@TODO see above @TODO>

    cluster_to_stats    #   can be ignored, it is used as a debugging tool


    USAGE:

    blacklist, cluster_to_stats = get_blacklisted_motif_clusters( e = 0.05, entropy_threshold = 10.0, species_threshold = 3.0, cluster_to_stats=None)

    """

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
    print '\tMean Entropy: '+str(H_mean)
    print '\tMean Number of Species: '+str(S_mean)
    print '\tNumber of Clusters: '+str(len(blacklist))+'         <-- no. of putative motifs\n'

    return H_mean, S_mean


#----------------------------
# clade grouping script
#----------------------------

def blacklist_to_clades( blacklist, cluster_to_stats, e ):

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


    **** Note that, consistent with vbname, e.g. MM10111_NGATTAAGCT_AGCTTAATCN_G, where '_G' can be:

        _G = motif is found in at least one gambiae complex species, 
        _A = motif is found in at least one anopheline species but not if that anopheline species is of the gambiae complex, 
        _M = motif is found in at least one mosquito species but not if that mosquito species is of the anophelines, 
        _D = motif is found in at least one dipteran species but not if that dipteran species is of the mosquitoes

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
    #   Here is where the logic is applied to group blacklisted motif clusters by Bob clades

    for c,d in blacklist:

        # GROUP THE @cluster INTO ONE OF THE CLADEs

        # DIPTERAN      # dipteran_clusters # contains at least one drosophila 
        if cluster_has_dipteran((c,d)): # store cluster if it contains at least one dipteran (drosophila)
            cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['bob_clade'] = 'D' # this allows us to @REFERENCE-BACK to which bobanian clade group this cluster belongs to, see: e.g.: CTRL+F: "def rename_blacklist(", note: D for dipteran, M for mosquito, A for anopheles, G for gambiae
            dipteran_clusters.append((c,d))

        # MOSQUITO      # mosquito_clusters # (contains no drosophila) AND (contains at least one culex OR one aedes)
        if ( cluster_has_mosquito((c,d)) and (not cluster_has_dipteran((c,d))) ): # store cluster if it contains at least one mosquito, but not if it has a dipteran
        # @TEST: any difference in "and not" vs. "and ( not ... "  ??
            cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['bob_clade'] = 'M' # see: CTRL+F: '@reference-back' for explanation
            mosquito_clusters.append((c,d))

        # ANOPHELES     # anopheles_clusters # (contains no mosquito AND no dipteran) AND (at least one anopheles (?) ) @BOB: wouldn't this just swallow up the whole of the gambiae complex...?..
        if cluster_has_anophelines((c,d)) and ( (not cluster_has_dipteran((c,d))) and (not cluster_has_mosquito((c,d)))):  # store cluster if it contains at least one anopheline, but not if it has a dipteran or a mosquito
            cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['bob_clade'] = 'A' # see: CTRL+F: '@reference-back' for explanation
            anopheles_clusters.append((c,d))

        # GAMBIAE       #..and the remaining 
        if cluster_has_gambiae_complex((c,d)) and ( (not cluster_has_dipteran((c,d))) and (not cluster_has_mosquito((c,d))) and (not cluster_has_anophelines((c,d)))): # store cluster if it contains at least one gambiae complex, but not if it has an anopheline or mosquito or dipteran
        
        #@TESTs_1: 
        #  - if the remainder == the output of this test, if yes, then great!
        #  - are there any overlaps between clusters in each of the groups? Hint: use set() 
            cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['bob_clade'] = 'G' # see: CTRL+F: '@reference-back' for explanation
            gambiae_clusters.append((c,d))

    print '\n_________________________________'
    print 'CLADE DISTRIBUTION OF BLACKLIST: '
    print '_________________________________\n'
    print '\tDipteran: '    +str(len(dipteran_clusters))
    print '\tMosquito: '    +str(len(mosquito_clusters))
    print '\tAnopheles: '   +str(len(anopheles_clusters))
    print '\tGambiae: '     +str(len(gambiae_clusters))
    print '\t_____________'
    print '\tTotal: '       +str(len(gambiae_clusters)+len(anopheles_clusters)+len(mosquito_clusters)+len(dipteran_clusters))
    print '\n'

    return dipteran_clusters, mosquito_clusters, anopheles_clusters, gambiae_clusters, cluster_to_stats


#----------------------------
# meme FBP outputs
#----------------------------

def make_meme_output(blacklist,outpath, cd_to_vbname, e=0.05):
    """

    Description:

    Take blacklisted clusters and generate a meme output file at a path specified by outpath arg

    Arguments:

    blacklist = blacklist, summary, clade_groups, cluster_to_stats = blacklist_then_summaryStats_from_shell(e,species_threshold,entropy_threshold)

    outpath = '../../data/stamp_data/out/dreme_100bp_e'+str(e)+'/blacklisted_motifs_meme_format.txt'

    """
    from time import sleep

    # blacklist --> this --> transfac and MEME matrices

    #-----------------------------
    # CLUSTER -to- FBP File path

    # Generates paths to FBP data for each cluster

    cluster_TO_fbp_path = {}

    for c,d in blacklist:

        name        = c
        distance    = d

        # If the c,d cluster motif was created from only a single (=1) dreme motif (i.e. whose name has ..._n001_..., it would have a different path to a motif cluster that was created from >1 dreme motifs. This is because STAMP will not create a *FBP.txt (PSSM) for single dreme motifs. This is because STAMP is software specialized for taking >1 motif and generating a consensus from them. So is the motif cluster generated from only a single dreme motif?

        # IF YES: 
        if 'n001' in name:

            datapath = '../../data/stamp_data/out/dreme_100bp_e'+str(e)+'/SWU_SSD/cluster_motifs_d'+str(distance)+'/fasta/'+str(name)+'.fasta' # e.g. 

        # IF NO:
        else:     
            datapath    = '../../data/stamp_data/out/dreme_100bp_e'+str(e)+'/SWU_SSD/cluster_motifs_d'+str(distance)+'/family_binding_matrices/'+str(name)+'FBP.txt' # e.g. ... c146_n004_TGGACGAKAFBP.txt # read tomtom inputdata format requirements

        cluster_TO_fbp_path[(c,d)] = datapath

    #-----------------------------
    # Blacklsited cluster -to- fbp TRANSFAC data

    cluster_TO_fbp = {}

    #print '\tTRANSFAC <Cluster Name> <Depth of Tree>' # headers for ...

    for c,d in blacklist:

        #print '\t'+str(c)+' '+str(d) # ... these values

        path = cluster_TO_fbp_path[(c,d)]

        cluster_TO_fbp[(c,d)] = {'transfac':[],'meme':[],'path':path}

        # IS there a *FBP.txt file available for the given c,d*?     *c,d = blacklisted cluster, where c is the cluster's name and d is the distance

        # YES:  then read the *FBP.txt file, incorporate the transfac motif PSSM (FBP), into a dict (see below)
        try:
            fi = open(path,'r')
            while True:
                l = fi.readline()
                if l == "":
                    break
                cluster_TO_fbp[(c,d)]['transfac'].append(l)
                #print l
            fi.close()
        # NO:  this means stamp has not generated a *FBP.txt (PSSM) for the motif: c,d ... this is because STAMP does not generate FBPs when there is only one input motif, e.g. all motifs whose names are ..._n001_... .
        except IOError:
            #sleep(5)
            print '\t\tWarning: No such file or directory: ../../data/stamp_data/out/dreme_100bp_e'+str(e)+'/SWU_SSD/cluster_motifs_d'+str(d)+'/family_binding_matrices/'+str(c)+'FBP.txt'

            print'\t\tWhat this means is this: cluster '+str(c)+'_d'+str(d)+' has a missing FBP that should have been generated by STAMP, and will thus count as missing data from downstream analyses...'

            cluster_TO_fbp[(c,d)]['transfac'].append('###MISSING###')

    #-----------------------------
    # Blacklsited cluster -to- fbp MEME data 

    for c,d in blacklist:

        fbp_transfac = cluster_TO_fbp[(c,d)]['transfac']

        # Is the (c,d) is missing it's FBP?
        if fbp_transfac == ['###MISSING###']:
        # YES: then we skip this motif cluster from downstream analyses 
            # cluster_TO_fbp[(c,d)]['meme']
            continue
        # NO:  then we do the following...
        else:
            #header  = 'MOTIF '+c+'_d_'+str(d)+'\n' # legacy motif names as generated by collapseMotifTree.py

            header = 'MOTIF '+cd_to_vbname[(c,d)]+'\n'   # convert a legacy name to vbname and store in "header" 

            rows = []

            for i,row in enumerate(fbp_transfac[1:-1]): # skips the header and tailing delimiter xx
                
                row_split = row.split('\t') # e.g. ['8', '0.0000', '0.0000', '1.0000', '0.0000', 'G\n']

                row_split_without_number_and_letter = row_split[1:-1] # e.g. ['0.0000', '0.0000', '1.0000', '0.0000']

                # A CHECK for if the FBP has it's PSSM in counts or decimals, and enforces decimals. e.g. [100,0,100,0] --> [0.5,0.0,0.5,0.0]
                row_as_floats = [float(x) for x in row_split_without_number_and_letter]
                row_as_total = sum(row_as_floats)
                if row_as_total > 1.1:
                    row_split_without_number_and_letter = [str(float(x)/float(row_as_total)) for x in row_as_floats]


                row_unsplit_without_number_and_letter = "\t".join(row_split_without_number_and_letter) # e.g. 0.3333\t0.3333\t0.3333\t0.0000


                rows.append(row_unsplit_without_number_and_letter)

            rows = "\n".join(rows)+"\n\n" # e.g. '0.3333\t0.3333\t0.3333\t0.0000\n0.0000\t0.0000\t1.0000\t0.0000\n1.0000\t0.0000\t0.0000\t0.0000\n1.0000\t0.0000\t0.0000\t0.0000\n0.1178\t0.0000\t0.0000\t0.8822\n0.3333\t0.6667\t0.0000\t0.0000\n0.0000\t0.1530\t0.0000\t0.8470\n1.0000\t0.0000\t0.0000\t0.0000\n0.0000\t0.0000\t1.0000\t0.0000\n\n'


            width=str(i+1)

            info    = 'letter-probability matrix: alength= 4 w= '+width+' nsites= 1000 E= 0.05\n' # nsites? E? 
            
            cluster_TO_fbp[(c,d)]['meme'] = header+info+rows
            # [header,info] + cluster_TO_fbp[(c,d)]['meme']

            # Reformat the fbp_transfac to make consistent w/ fbp_meme
            fbp_transfac_str = "\n".join(fbp_transfac)
            cluster_TO_fbp[(c,d)]['transfac'] = fbp_transfac_str


    #-----------------------------
    # Write to file

    print '\t\t'+outpath

    fo = open(outpath,'w')
    fo.write("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.303 C 0.183 G 0.209 T 0.306 \n\n")
    for c,d in cluster_TO_fbp.keys():
        if cluster_TO_fbp[(c,d)]['meme']: # If statement here, because some of the (c,d)s have missing (empty) values such as [], and this will crash the fo.write(), so we skip the (c,d)s with empty values.
            fo.write(cluster_TO_fbp[(c,d)]['meme'])
    fo.close()

def generate_summary_statistics_file(species_threshold, entropy_threshold, outDirPath, summary, blacklist, clade_groups, cd_to_vbname, cluster_to_stats, e = 0.05):

    """

    Description:

    Generates a file with summary statistics of a given blacklist, in a similar fashion to: blacklist_to_clades() and print_blacklist_summary(). It also provides the list of motif clusters of the blacklist and separate lists according to Bob-clade-groupings.

    Arguments:

    ...

    Tags:

    @stats


    """

    def fraction_of_motifs_matching_to_jaspar( match_significance, blacklist ):

        """

        Returns the fraction of motifs in the blacklist with a match in JASPAR for a given p-value threshold that defines what constitues a match. If the p-value of a match is 0.05 then the fraction is calculated by the number of motifs in the blacklist with at least one JASPAR match whose p-value is less than the threshold.

        """

        n_jaspar_matched_motifs = 0 

        for c,d in blacklist:

            if cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['cluster']['data']['jaspar']['stats']['<=p'][match_significance]>0:
                n_jaspar_matched_motifs += 1

        try:
            fraction_of_motifs_matched = float(n_jaspar_matched_motifs)/len(blacklist)
        except ZeroDivisionError:
            fraction_of_motifs_matched = 0 # ANDY: sometimes len(blacklist) = 0, for now I deal with this by making JASPAR statistic = 0

        return fraction_of_motifs_matched


    ###################################################################
    # Store JASPAR matches per motif of a collpased motif cluster c,d
    ###################################################################

    for c,d in blacklist: # e.g. c047_n007_CATTGGGCG, 0.29799999999999999

        # Opportunistic corrections to cluster_to_stats
        cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['cluster']['path']['jaspar'] = cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['cluster']['path']['jaspar'].replace('paris','pairs')
        cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['cluster']['data']['jaspar'] = {'matches':{},'stats':{}} # represents matches made per motif of the collapsed motif cluster


        cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['cluster']['data']['jaspar']['stats']['n_total'] = cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['motif']['n']
        cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['cluster']['data']['jaspar']['stats']['n_queries'] = 0


        # access *_match_pairs.txt file
        path    = '../'+cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['cluster']['path']['jaspar'] # e.g. '../../data/stamp_data/out/dreme_100bp_e0.05/SWU_SSD/cluster_motifs_d0.298/family_binding_matrices/c047_n007_CATTGGGCG_match_pairs.txt'
        fi      = open(path,'r')

        # initiate
        matches     = []
        query       = []
        

        while True:

            line = fi.readline() # e.g. Snail   5.7055e-01  CCGCCACTG   ---CACCTG

            if line == '':
                break

            if '>' in line: # e.g. >   AAEG_GTKTAGTGA
                # STORE QUERY:[MATCH1,MATCH2,...]
                cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['cluster']['data']['jaspar']['stats']['n_queries'] += 1 # update counting of queries with matches
                if matches:
                    cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['cluster']['data']['jaspar']['matches'][query] = matches
                # RESET QUERY:MATCHES
                query   = line.rstrip().split('\t')[1] # e.g. AAEG_GTKTAGTGA
                matches = []
            else:
                split_line      = line.rstrip().split('\t')
                match_name      = split_line[0] # e.g. Myf
                match_pvalue    = float(split_line[1]) # e.g. 2.9990e-01
                match_qseq      = split_line[2] # e.g. CCGCCACTG---
                match_mseq      = split_line[3] # e.g. MRGCARCTGCTG
                matches.append([match_name,match_pvalue,match_qseq,match_mseq]) # e.g. ['Myf',2.9990e-01,'CCGCCACTG---','MRGCARCTGCTG']

        fi.close()



    #########################################################################################################################
    # Count the number of query motifs with a p-value <= p-value threshold to each (c,d), for various p-value thresholds
    ##########################################################################################################################


    p_thresh_vec = [0.05,0.01,0.001] # significance thresholds: 0.05, 0.01, 0.001

    for c,d, in blacklist:

        cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['cluster']['data']['jaspar']['stats']['<=p'] = {}

        for p_thresh in p_thresh_vec:

            p_thresh_n_match    = 0

            for query in cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['cluster']['data']['jaspar']['matches'].keys():

                matches             = cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['cluster']['data']['jaspar']['matches'][query]

                p_thresh_matches    = [] # number of queries with at least one match

                for match in matches:
                    
                    pvalue  = match[1]

                    if pvalue <= p_thresh:
                        
                        p_thresh_matches.append(match)

                if p_thresh_matches: # if there are any matches with p-values <= p-value threshold, then +1 

                    p_thresh_n_match += 1

            cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['cluster']['data']['jaspar']['stats']['<=p'][p_thresh] = p_thresh_n_match


    #
    # JASPAR STATISTIC for significance thresholds: 0.05, 0.01, 0.001  
    # 

    jaspar_statistics=[]

    for match_significance in p_thresh_vec:

        jaspar_statistic = fraction_of_motifs_matching_to_jaspar( match_significance, blacklist )

        jaspar_statistics.append(str(jaspar_statistic)+' (matches p-value<='+str(match_significance)+')')





    fo = open(outDirPath+'full_summary_statistics.txt','w')


    # SUMMARY STATISTICS 

    H_mean = summary[0]
    S_mean = summary[1]

    fo.write('_________________________________\n')
    fo.write('SUMMARY STATISTICS OF BLACKLIST: \n')
    fo.write('_________________________________\n')
    fo.write('\tMean Entropy: '+str(H_mean)+'\n')
    fo.write('\tMean Number of Species: '+str(S_mean)+'\n')
    fo.write('\tNumber of Clusters: '+str(len(blacklist))+'\n')
    fo.write('\n\n')

    fo.write('_________________________________\n')
    fo.write('JASPAR STATISTICS (the fraction of motifs in blacklist with matches to JASPAR): \n')
    fo.write('_________________________________\n')
    for jaspar_statistic in jaspar_statistics:
        fo.write('\t'+jaspar_statistic+'\n')
    fo.write('\n\n')

    dipteran_clusters   = clade_groups[0]
    mosquito_clusters   = clade_groups[1]
    anopheles_clusters  = clade_groups[2]
    gambiae_clusters    = clade_groups[3]

    fo.write('_________________________________\n')
    fo.write('CLADE DISTRIBUTION OF BLACKLIST: \n')
    fo.write('_________________________________\n')
    fo.write('\tDipteran: '    +str(len(dipteran_clusters))+'\n')
    fo.write('\tMosquito: '    +str(len(mosquito_clusters))+'\n')
    fo.write('\tAnopheles: '   +str(len(anopheles_clusters))+'\n')
    fo.write('\tGambiae: '     +str(len(gambiae_clusters))+'\n')
    fo.write('\t_____________\n')
    fo.write('\tTotal: '       +str(len(gambiae_clusters)+len(anopheles_clusters)+len(mosquito_clusters)+len(dipteran_clusters))+'\n')
    fo.write('\n\n\n')


    # NAMES OF CLUSTERS LISTED by the following categories: whole blacklist, dipteran, mosquito, anopheles, gambiae

        # FULL BLACKLIST

    fo.write('Names of putative motifs (collapsed clusters of dreme motifs) grouped by:\n\n')

    fo.write('all:\n')
    fo.write('<VB Name> \t <Old Name> \t <Group> \t <Entropy> \t <nspecies> \t <nparalogues> \t <avglength>\n') # headers to indicate format of motif names

    for c,d in blacklist:

        ###########################
        # per cluster STATISTICS 
        ###########################

        vbname      = cd_to_vbname[(c,d)] # e.g. MM10001_CAKTGGCGG_CCGCCAMTG_M
        group       = cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['bob_clade']
        entropy     = cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['cluster']['H']

        if str(entropy) == '-0': # for some reason, some entropies are: '-0' ... so we correct it
            entropy = 0.0

        nspecies    = cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['species']['unique']['n_unique']
        nparalogues = cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['species']['unique']['n_paralogue']
        avglength   = cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['motif']['avg_length']

        fo.write(vbname+'\t'+c+'_d'+str(d)+'\t'+str(group)+'\t'+str(entropy)+'\t'+str(nspecies)+'\t'+str(nparalogues)+'\t'+str(avglength)+'\n') # e.g. MM10001_CAKTGGCGG_CCGCCAMTG_M   c047_n007_CATTGGGCG_d0.298

    fo.write('\n\n')


        # PER CLADE 

    fo.write('Blacklisted motifs grouped by clade:\n')

    i_TO_clade = {0:'dipteran',1:'mosquito',2:'anopheles',3:'gambiae'}

    for i,blacklist_clade in enumerate(clade_groups):

        clade = i_TO_clade[i] # e.g. 'dipteran'

        fo.write('\n'+clade+':\n') 
        fo.write('<Motif_name> \t <Distance_in_motif_tree>\n')

        for c,d in blacklist_clade:

            # is there is a translation available for (c,d) to vb name?
            vbname = cd_to_vbname[(c,d)]

            fo.write(vbname+'\t'+str(c)+'_'+str(d)+'\n') # e.g. i2 = MM10001_CAKTGGCGG_CCGCCAMTG_M  c047_n007_CATTGGGCG_d0.298

        fo.write('\n\n')

    fo.close()


def rename_blacklist( blacklist, cluster_to_stats, e = 0.05, version = 1 ):

    print 'Re-naming blacklisted motif clusters (c,d) to vb-name format...'

    cd_to_vbname = {}   # @dict_1 that converts (c,d) to vbname..

    #   ..general format: MM<version><id, e.g. 0001>_<STAMP ambiguity seq>_<reverse complement>..
    #   ..e.g. MM10001_CTGATGGNC_GNCCATCAG 

    version = int(version)             # @WARN: this version number will need to change in future motif versions


    for n,(c,d) in enumerate(blacklist):   # e.g. (c,d) = ('c397_n002_CTGATGGSC','0.056')

        # grab the collapsed motif's STAMP (or .fasta) motif sequence

        fbp_path = cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['cluster']['path']['fbp'] # e.g. '../data/stamp_data/out/dreme_100bp_e0.05/SWU_SSD/cluster_motifs_d0.056/family_binding_matrices/c397_n002_CTGATGGSCFBP.txt'
        fbp_fi  = open('../'+fbp_path,'r') 

        # ^ e.g.                              \/----letters of seq-----\/
        # .._n000(>1)_.. :                          .._n0001_.. :
        # DE  FBP                                 # DE CQUI_ATATGGAGA
        # 0   0.1250  0.6250  0.1250  0.1250  c   # 0   198 0   0   0   A
        # 1   0.0000  0.0000  0.2178  0.7822  T   # 1   0   0   0   198 T
        # 2   0.0000  0.0000  1.0000  0.0000  G   # 2   198 0   0   0   A
        # 3   0.7513  0.0000  0.2487  0.0000  A   # 3   0   0   0   198 T
        # 4   0.0000  0.0000  0.0000  1.0000  T   # 4   0   0   198 0   G
        # 5   0.0000  0.0000  1.0000  0.0000  G   # 5   0   0   198 0   G
        # 6   0.0000  0.0000  1.0000  0.0000  G   # 6   198 0   0   0   A
        # 7   0.1250  0.4198  0.3302  0.1250  n   # 7   0   0   198 0   G
        # 8   0.1250  0.6250  0.1250  0.1250  c   # 8   198 0   0   0   A
        # XX                                      # XX

        seq = []                           #  ^-----e.g.----------------^        
        while True:
            l = fbp_fi.readline()
            if l == "XX\n": break
            if l == 'DE\tFBP\n': continue
            s = l[-2] # e.g. n
            seq.append(s) # e.g. ['c', 'T', 'G', 'A', 'T', 'G', 'G', 'n', 'c']
        fbp_fi.close()

        # find the complement of seq,   
        #   e.g. if seq = ['c', 'T', 'G', 'A', 'T', 'G', 'G', 'n', 'c'] then rc = ['G', 'N', 'C', 'C', 'A', 'T', 'C', 'A', 'G']

        iupac_to_complement = {'A':'T','C':'G','G':'C','T':'A','M':'K','R':'Y','W':'W','S':'S','Y':'R','K':'M','V':'B','H':'D','D':'H','B':'V','N':'N'}
        complement = [iupac_to_complement[i.upper()] for i in seq]
        rc = ''.join(complement[::-1]) # e.g. CCGCCAMTG

        # format to strings
        seq = ''.join(seq).upper()
        rc = ''.join(rc).upper()

        # Grab the 'bob_clade' label, e.g. 'D' for dipteran, 'M' for mosquito, etc.
        bob_clade_symbol = cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['bob_clade']

        # vbname = 'MM'+str(version)+str(n+1).zfill(4)+'_'+seq+'_'+rc+'_'+str(bob_clade_symbol) # e.g. MM10001_CTGATGGNC_GNCCATCAG_D

        vbname = 'MM'+str(version)+str(n+1).zfill(4)+'_'+seq+'_'+rc # e.g. MM10001_CTGATGGNC_GNCCATCAG ^ @pad 0s to the left to have 4 digits, e.g. 0001 and 0111, etc. so that we have a max of 9999 motifs we can use, keeping motif names at controlled numbering width.
        #                           ^ n+1 to ensure indexing starts from 1 (not 0)

        # incorporate into @dict_1
        cd_to_vbname[(c,d)] = vbname


    return cd_to_vbname


def blacklist_then_summaryStats_S_and_H_combos( S_from, S_to, S_step , H_from, H_to, H_step , e = 0.05):

    """

    Description:

        runs blacklist_then_summaryStats_from_shell() on combinations of -nspecies, -entropy parameters

        Output is found in: /home/ab108/0VB/2kb/data/stamp_data/out/dreme_100bp_e0.05/SWU_SSD/progressively_collapsed_motifs/

    Arguments:

        ...@todo...

    Usage:

        blacklist_then_summaryStats_S_and_H_combos( S_from = 1, S_to = 21, S_step = 1, H_from = 1, H_to = 30, H_step = 1, e = 0.05)

    Tags:

        @combo

    """

    import numpy as np
    from itertools import product

    S_vec = np.arange(float(S_from),float(S_to),float(S_step))  # vector of -nspecies thresholds
    H_vec = np.arange(float(H_from),float(H_to),float(H_step))  # vector of -entropy thresholds
    S_and_H_vector = list(product(S_vec,H_vec))  # combinations of -nspecies, -entropy thresholds

    cluster_to_stats = load_motif_cluster_stats_dict() # load data


    print 'Threshold Parameters: '

    for S,H in S_and_H_vector:

        print 'Number of species >= : '+str(S)+'  Entropy <= : '+str(H)+'\n\n'
     
        blacklist, summary, clade_groups, cluster_to_stats, cd_to_vbname = blacklist_then_summaryStats_from_shell( e = 0.05, species_threshold = S, entropy_threshold = H, cluster_to_stats=cluster_to_stats )







#----------------------------
# MAIN THING
#----------------------------

def blacklist_then_summaryStats_from_shell( e = 0.05, species_threshold = 3, entropy_threshold = 5,cluster_to_stats=None ):
    """ 

    Description:\n\n

        A wrapper for obtaining the blacklist at given species number threshold, S, and entropy threshold, H, (CTRL+F: blacklist_criteria_check) and e value (CTRL+F: @e-value)\n\n

    Arguments:\n\n 

        e = 0.05                      # CTRL+F: @e-value\n\n

        species_threshold = 3.0       # CTRL+F: blacklist_criteria_check\n\n

        entropy_threshold = 1.0       # CTRL+F: blacklist_criteria_check\n\n

        Example:\n\n

        blacklist_then_summaryStats_from_shell( e = 0.05, species_threshold=3.0, entropy_threshold=1.0)\n\n

    Usage:

        blacklist, summary, clade_groups, cluster_to_stats, cd_to_vbname = blacklist_then_summaryStats_from_shell( e = 0.05, species_threshold = 3, entropy_threshold = 5,cluster_to_stats=None )

    alt Usage:

        blacklist, summary, clade_groups, cluster_to_stats, cd_to_vbname = blacklist_then_summaryStats_from_shell( e = 0.05, species_threshold = 3, entropy_threshold = 5,cluster_to_stats=cluster_to_stats ) # if cluster_to_stats.p is already loaded..

    Tags:

        @shell


    """
    import os
    import shutil

    # GENERATE @blacklist

    if cluster_to_stats == None: # Skips the cluster_to_stats loading if it is already assigned
        blacklist, cluster_to_stats = get_blacklisted_motif_clusters( e, entropy_threshold, species_threshold,  cluster_to_stats=None) 
    else:
        blacklist, cluster_to_stats = get_blacklisted_motif_clusters( e, entropy_threshold, species_threshold,  cluster_to_stats=cluster_to_stats)


    # GENERATE cd_to_vbname translations


    print '\nGenerating statistics...\n'

    # GLOBAL STATISTICS
    summary = print_blacklist_summary( blacklist, cluster_to_stats, e)# (H_mean, S_mean)

    # BOB CLADE GROUPINGS
    
    clade_group_data = blacklist_to_clades( blacklist, cluster_to_stats, e ) # dipteran_clusters, mosquito_clusters, anopheles_clusters, gambiae_clusters, cluster_to_stats

    clade_groups = clade_group_data[:-1]# dipteran_clusters, mosquito_clusters, anopheles_clusters, gambiae_clusters
    
    clade_groups_with_all = list(clade_groups) + [blacklist]  # append the whole blacklist so we can use it as "all"

    cluster_to_stats = clade_group_data[-1] # cluster_to_stats // note: this has been modified now to allow referencing to the clade groupings explicitly via cluster_to_stats, e.g. cluster_to_stats['e'+str(e)+'_d'+str(d)][c]['bob_clade'] = 'D'     }- where 'D' means dipteran

    cd_to_vbname = rename_blacklist( blacklist, cluster_to_stats, e=e, version = 1 ) # dictionary: converts a blacklisted motif cluster (c,d) into a vectorbase name e.g. of usage: vbname = cd_to_vbname[(c,d)] // e.g. of vbname = 'MM10119_CGATGCGAT_ATCGCATCG' // general format: MM<version><id, e.g. 0001>_<STAMP ambiguity seq>_<reverse complement>..

    i_TO_clade = {0:'dipteran',1:'mosquito',2:'anopheles',3:'gambiae',4:'all'}
    
    print 'Generating MEME output data per clade, apply TOMTOM on these...\n'

    # GENERATE MEME DATA FILE SYSTEM

    if not os.path.exists('../../data/stamp_data/out/dreme_100bp_e'+str(e)+'/SWU_SSD/progressively_collapsed_motifs/'):
        os.makedirs('0VB/2kb/data/stamp_data/out/dreme_100bp_e'+str(e)+'/SWU_SSD/progressively_collapsed_motifs/') # generate parent directory that stores data for all --nspecies, --entropy combinations

    outdir = '../../data/stamp_data/out/dreme_100bp_e'+str(e)+'/SWU_SSD/progressively_collapsed_motifs/nspecies_'+str(species_threshold)+'_entropy_'+str(entropy_threshold)+'/'
    
    if os.path.exists(outdir): # remove an output directory if it exists, i.e. overwrite them, e.g. of an outdir is /nspecies_3.0_entropy_5.0/
        shutil.rmtree(outdir)
    os.makedirs(outdir)

    generate_summary_statistics_file(species_threshold, entropy_threshold, outdir, summary, blacklist, clade_groups, cd_to_vbname, cluster_to_stats, e=e) # Generates a file with verbose summary statistics

    for i,blacklist_clade in enumerate(clade_groups_with_all):
        clade = i_TO_clade[i]
        print '\t'+clade
        
        # OVERWRITE directory to house the MEME data, e.g. /SWU_SSD/nspecies_3_entropy_10/

        # TODO: <some function that generates the list of blacklisted motif cluster names + four other lists of cluster names: one per bobainian clade group>

        outpath = outdir+clade+'_blacklisted_meme_format.txt'

        make_meme_output(blacklist_clade,outpath, cd_to_vbname, e) # generate the meme FBP data

    return blacklist, summary, clade_groups, cluster_to_stats, cd_to_vbname



#---------------------------------------------------------
# Scraps




#---------------------------------------------------------




import sys,getopt


#@WARN: this script can be run from termnial, but NEEDS two arguments to be fed to it: (i) from terminal: --nspecies 3.0 --entropy 10.0, (ii) from Python: blacklist, summary, clade_groups, cluster_to_stats = blacklist_then_summaryStats_from_shell( e = 0.05, species_threshold = 3.0, entropy_threshold = 10.0 )

# Argument Handling
def main(argv):

    # USER MADE ERROR? 

    # ... if no
    try:

        opts, args = getopt.getopt(argv,"hs:H",["nspecies=","entropy="])

    # ... if yes (e.g. mispelled args, non-existant args etc)
    except getopt.GetoptError:
        
        print 'collapseMotifTree_progressiveMode.py --nspecies <species threshold,e.g.3> --entropy <entropy threshold,e.g.10.5>'

        sys.exit(2)



    # PROCESS THE @ARGS

    for opt, arg in opts:
        
        if opt == '-h': # help file 
        
            print 'collapseMotifTree_progressiveMode.py -s <species threshold> -H <entropy threshold>\n\nDescription:\n\nA wrapper for obtaining the blacklist at given species number threshold, S, and entropy threshold, H, (CTRL+F: blacklist_criteria_check) and e value (CTRL+F: @e-value)\n\nArguments:\n\n e = 0.05                      # CTRL+F: @e-value\n\nspecies_threshold = 3.0       # CTRL+F: blacklist_criteria_check\n\nentropy_threshold = 1.0       # CTRL+F: blacklist_criteria_check\n\nExample:\n\nblacklist_then_summaryStats_from_shell( e = 0.05, species_threshold=3.0, entropy_threshold=1.0)\n\n'

            sys.exit()

        elif opt in ("-s","--nspecies"):    # long and short args, see

            s = arg # don't forget to correct the type!!

        elif opt in ("-H", "--entropy"):

            H = arg

    # CMD or PYTHON EXECUTE CONTROL: by default try run as if on cmd line, with input and output args ... else if it errors it may be running from within python, so then run default args: 

    # ... script was run from command line
    try:

        # SCRIPT
        blacklist_then_summaryStats_from_shell( e = 0.05, species_threshold = float(s), entropy_threshold = float(H) ) # <-- main script, @REPLICATE exactly as below
        #                            ^@args here processed from above

    # ... script was run within Python, e.g. %paste in IPython or imported
    except UnboundLocalError:
        
        s = raw_input('species threshold? e.g. s = 3.0...')
        
        H = raw_input('entropy threshold? e = 10.0...')

        blacklist, summary, clade_groups, cluster_to_stats =blacklist_then_summaryStats_from_shell( e = 0.05, species_threshold = float(s), entropy_threshold = float(H) ) # <-- main script, CTRL+F "@replicate" exactly as above

if __name__ == "__main__":
    main(sys.argv[1:])


