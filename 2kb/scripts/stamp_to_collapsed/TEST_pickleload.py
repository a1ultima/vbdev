# NOTES:
# Some technical words may not make sense to you. Such words may have intuitive descriptions available in here. Words with intuitive descriptions available will be in the format: @<word-in-lowercase>. Simply CTRL+F "@<word-in-uppercase>" to jump to the intuitive descriptions.


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
                                                     # WARN: should test this, got it from StackOverflow: http://stackoverflow.com/a/3170067/3011648
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

    for c in blacklist:

        # GROUP THE @cluster INTO ONE OF THE CLADEs

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
    print '\t_____________'
    print '\tTotal: '       +str(len(gambiae_clusters)+len(anopheles_clusters)+len(mosquito_clusters)+len(dipteran_clusters))
    print '\n'

    return dipteran_clusters, mosquito_clusters, anopheles_clusters, gambiae_clusters


#----------------------------
# meme FBP outputs
#----------------------------

def make_meme_output(blacklist,outpath,e=0.05):
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
            header  = 'MOTIF '+c+'_d_'+str(d)+'\n'

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

def generate_summary_statistics_file(species_threshold, entropy_threshold, outDirPath, summary, blacklist, clade_groups):

    """

    Description:

    Generates a file with summary statistics of a given blacklist, in a similar fashion to: blacklist_to_clades() and print_blacklist_summary(). It also provides the list of motif clusters of the blacklist and separate lists according to Bob-clade-groupings.

    Arguments:

    ...

    """

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
    fo.write('<Motif_name> \t <Distance_in_motif_tree>\n') # headers to indicate format of motif names

    for c in blacklist:

        fo.write(c[0]+'\t'+str(c[1])+'\n')

    fo.write('\n\n')


        # PER CLADE 

    fo.write('Blacklisted motifs grouped by clade:\n')

    i_TO_clade = {0:'dipteran',1:'mosquito',2:'anopheles',3:'gambiae'}

    for i,blacklist_clade in enumerate(clade_groups):

        clade = i_TO_clade[i]

        fo.write('\n'+clade+':\n')
        fo.write('<Motif_name> \t <Distance_in_motif_tree>\n')

        for i2 in blacklist_clade:
            fo.write(str(i2[0])+'\t'+str(i2[1])+'\n')

        fo.write('\n\n')

    fo.close()


#----------------------------
# MAIN THING
#----------------------------

def blacklist_then_summaryStats_from_shell( e = 0.05, species_threshold = 3.0, entropy_threshold = 1.0,cluster_to_stats=None ):
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
    import os
    import shutil

    # GENERATE @blacklist

    if cluster_to_stats == None: # Skips the cluster_to_stats loading if it is already assigned
        blacklist, cluster_to_stats = get_blacklisted_motif_clusters( e, entropy_threshold, species_threshold,  cluster_to_stats=None) 
    else:
        blacklist, cluster_to_stats = get_blacklisted_motif_clusters( e, entropy_threshold, species_threshold,  cluster_to_stats=cluster_to_stats)


    print '\nGenerating statistics...\n'

    # GLOBAL STATISTICS
    summary = print_blacklist_summary( blacklist, cluster_to_stats, e)# (H_mean, S_mean)

    time.sleep(2)

    # BOB CLADE GROUPINGS
    clade_groups = blacklist_to_clades( blacklist, cluster_to_stats, e ) # dipteran_clusters, mosquito_clusters, anopheles_clusters, gambiae_clusters
    i_TO_clade = {0:'dipteran',1:'mosquito',2:'anopheles',3:'gambiae'}
    
    time.sleep(2)


    
    print 'Generating MEME output data per clade, apply TOMTOM on these...\n'

    # GENERATE MEME DATA FILE SYSTEM

    if not os.path.exists('../../data/stamp_data/out/dreme_100bp_e'+str(e)+'/SWU_SSD/progressively_collapsed_motifs/'):
        os.makedirs('../../data/stamp_data/out/dreme_100bp_e'+str(e)+'/SWU_SSD/progressively_collapsed_motifs/') # generate parent directory that stores data for all --nspecies, --entropy combinations

    outdir = '../../data/stamp_data/out/dreme_100bp_e'+str(e)+'/SWU_SSD/progressively_collapsed_motifs/nspecies_'+str(species_threshold)+'_entropy_'+str(entropy_threshold)+'/'
    
    if os.path.exists(outdir): # remove an output directory if it exists, i.e. overwrite them, e.g. of an outdir is /nspecies_3.0_entropy_5.0/
        shutil.rmtree(outdir)
    os.makedirs(outdir)

    generate_summary_statistics_file(species_threshold, entropy_threshold, outdir, summary, blacklist, clade_groups) # Generates a file with verbose summary statistics

    for i,blacklist_clade in enumerate(clade_groups):
        clade = i_TO_clade[i]
        print '\t'+clade
        
        # OVERWRITE directory to house the MEME data, e.g. /SWU_SSD/nspecies_3_entropy_10/

        # TODO: <some function that generates the list of blacklisted motif cluster names + four other lists of cluster names: one per bobainian clade group>

        outpath = outdir+clade+'_blacklisted_meme_format.txt'

        make_meme_output(blacklist_clade,outpath,e) # generate the meme FBP data

    return blacklist, summary, clade_groups, cluster_to_stats


#---------------------------------------------------------
# Scraps




#---------------------------------------------------------





import numpy as np
from itertools import product
import time


cluster_to_stats = load_motif_cluster_stats_dict() # load data



# blacklist, summary, clade_groups, cluster_to_stats = blacklist_then_summaryStats_from_shell( e = 0.05, species_threshold = 1.0, entropy_threshold = 5.0, cluster_to_stats=cluster_to_stats )






# combinations of -nspecies, -entropy

e = 0.05 # 
S_vec = np.arange(1,21,1)  # vector of -nspecies thresholds
H_vec = np.arange(3,30,1)  # vector of -entropy thresholds
S_and_H_vector = list(product(S_vec,H_vec))  # combinations of -nspecies, -entropy thresholds


print 'Threshold Parameters: '

for S,H in S_and_H_vector:

    print 'Number of species >= : '+str(S)+'  Entropy <= : '+str(H)+'\n\n'
 
    blacklist, summary, clade_groups, cluster_to_stats = blacklist_then_summaryStats_from_shell( e = 0.05, species_threshold = S, entropy_threshold = H, cluster_to_stats=cluster_to_stats )







