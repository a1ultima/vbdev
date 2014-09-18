import os
from ete2 import Tree
import ete2
import numpy as np 
import matplotlib.pyplot as plt
import datetime
import subprocess
import math


#README
#   TODO: still needs to be transformed into a robust and proper function, with commandline args

"""

Description:

...

Run the following scripts in order, in the way shown:

1. collapseMotifTree.py # set the parameters at @PARAMETERS or @WARN
2. motifStatistics.py
3. python collapseMotifTree_progressiveMode.py --nspecies 3.0 --entropy 10.0

Notes:

CTRL+F: "@WARN" for the warnings

"""



# SUB-FUNCTIONS:
def mean(array):
    return sum(array)/float(len(array))

def cache_distances(tree):
    ''' precalculate distances of all nodes to the root''' 
    node2rootdist = {tree:0}
    for node in tree.iter_descendants('preorder'):
        node2rootdist[node] = node.dist + node2rootdist[node.up]
    return node2rootdist

def collapse(tree, min_dist):
    # cache the tip content of each node to reduce the number of times the tree is traversed
    node2tips       = tree.get_cached_content()
    root_distance   = cache_distances(tree)
    for node in tree.get_descendants('preorder'):
        if not node.is_leaf():
            avg_distance_to_tips = mean([root_distance[tip]-root_distance[node] for tip in node2tips[node]])
            if avg_distance_to_tips < min_dist:
                # do whatever, ete support node annotation, deletion, labeling, etc.
                # rename
                node.name += ' COLLAPSED avg_d:%g {%s}' %(avg_distance_to_tips,','.join([tip.name for tip in node2tips[node]]))
                # label
                node.add_features(collapsed=True)
                # set drawing attribute so they look collapsed when displayed with tree.show()
                node.img_style['draw_descendants'] = False
                # etc...

# MAIN FUNCTION
def exploreCutoffs( dataPath, filename_in, d_from = 0, d_to = 2, d_step = 0.05, save_d = 0.5 ):
    """
    Description:

    Given a Newick format tree, collapses clades* whose lengths to tips are on average < d, where d is here incremented from e.g. 0 to 1 in 0.05 increments.

    *in the case of this project, clades represent clusters of motifs that are collapsed together and their PSSMs merged. These will be later referred to as @motif clusters

    Arguments:

    dataPath    = '../../data/stamp_data/SWU_SSD/'  # directory of motif tree
    filename_in = 'out.tree'                        # Newick format motif tree
    d_from      = 0.00                              # distance cut-off start
    d_to        = 1.00                              # distance cut-off end
    d_step      = 0.05                              # increment of distance cut-off
    save_d      = 0.5                               # save the clustered motif names if the loop cuts at a d within save_d
    """

    # TEST  :   catch annoying errors where the distance cutting range (d_from, d_to) does not include the save_d !!!!
    #if not str(save_d) in str(np.arange(d_from,d_to,d_step)):

    if not str(save_d) in str(np.arange(d_from,d_to,d_step)):

        raise Exception(" save_d does not lie in the distance-cut-off range!! this is asking for trouble... you won't get your tree saved... look: "+"Tree Range: ")

    print('initiating tree cutting...')

# Statistics and data vars
    clusterTree_arr         = []
    clusterSize_dist_arr    = []
    clusterSize_mean_arr    = []
    clusterSize_median_arr  = []
    clusterSize_std_arr     = []
    distances               = np.arange(d_from,d_to,d_step)# generate distance vector
    t                       = Tree(dataPath+filename_in)    # Load a tree structure from a newick file.

    # MAIN LOOP     :   for each distance cut-off make a collapse
    print 'Trees are collapsing...'

    for count,d in enumerate(distances):     # explore distance threshold params
        
        if count%1 == 0:
            print( 'Tree: '+str(count)+'  Distance: '+str(save_d) )

        t = Tree(dataPath+filename_in)           # load tree

        collapse(t, d)                           # collapse /u threshold d
        
        for n in t.search_nodes(collapsed=True): # re-format tree via collapses
            for ch in n.get_children():
                ch.detach()

        # node distances
        dist_dict = cache_distances(t)           # recorbranch lengths of newly collapsed tree

        # unite nested list of collapsed motifs
        clusters_multiple   = [i.name.split('{')[1][:-1].split(',') for i in dist_dict.keys() if "COLLAPSED" in i.name] # list of clusters, where each cluster is a list of motifs)
        clusters_single     = [[i.name] for i in dist_dict.keys() if not 'NoName' in i.name]
        clusters_all        = clusters_single+clusters_multiple

        # SAVE CLUSTER NAMES, HEADERS_FOR_FASTA  :   saves on a specified collapsing distance cut-off, save_d
        #print(str(save_d),str(d))

        if str(save_d) == str(d): # @SAVED

            save_dir = dataPath+'cluster_motifs_d'+str(save_d)+'/headers_for_fasta/' # path to store files

            # Remove headers_for_fasta folder
            if os.path.exists(save_dir):
                print('clustered tree folder already exists, overwriting...')
                import shutil
                shutil.rmtree(save_dir)
                os.makedirs(save_dir)
            else:
                os.makedirs(save_dir)
            
            for i,cluster in enumerate(clusters_all):
                try: # works for clusters_multiple (nested list)
                    filename_cluster = 'c'+str(i).zfill(3)+'_n'+str(len(cluster)).zfill(3)+'_'+cluster[0][0].split('_')[1]+'.txt'
                except IndexError: # works for clusters_single (flat list)
                    filename_cluster = 'c'+str(i).zfill(3)+'_n'+str(len(cluster)).zfill(3)+'_'+cluster[0].split('_')[1]+'.txt'
                file_out_cluster = open(save_dir+filename_cluster,'w')
                for motif in cluster:
                    file_out_cluster.write(motif+'\n')

    # ranked distribution of cluster sizes, largest-to-smallest
        clusterTree         = t
        clusterSize_dist    = np.array(sorted([float(len(i)) for i in clusters_all],reverse=True))
        clusterSize_mean    = np.mean(clusterSize_dist) 
        clusterSize_median  = np.median(clusterSize_dist)
        clusterSize_std     = np.std(clusterSize_dist)

    # store
        clusterTree_arr.append(clusterTree)
        clusterSize_dist_arr.append(clusterSize_dist) 
        clusterSize_mean_arr.append(clusterSize_mean)
        clusterSize_median_arr.append(clusterSize_median)
        clusterSize_std_arr.append(clusterSize_std)

# write collapsed newicks
    for j,i in enumerate(clusterTree_arr):
        if not os.path.exists(dataPath+'collapsed_newick/'):
            os.makedirs(dataPath+'collapsed_newick/')
        file_out = open(dataPath+'collapsed_newick/d'+str(distances[j])+'.tre','w')
        file_out.write(i.write())
        file_out.close()
        # print 'd: '+str(distances[j])+'  '+str(i.write())

# UNCOMMENT FOR SUMMARY STATS FOR THE d_CUTOFFS
#         # print '   '
# # Vectors to plot
#     x   = distances
#     y   = clusterSize_mean_arr
#     z   = clusterSize_median_arr
#     yerr= [0]*len(distances)
#     xerr= clusterSize_std_arr

    # # diagram
    # plt.errorbar(x, y, xerr, yerr, capsize=0, ls='none', color='gray', elinewidth=1)
    # plt.annotate('Std deviation bars', (0.55, 58), xytext=(-10, 60), va='bottom',textcoords='offset points',arrowprops=dict(facecolor='gray', shrink=0.05),color='Gray')

    # # Text to annotate on plot
    # ydescrip = [str(int(round(i))) for i in y] # means 
    # zdescrip = [str(int(round(i))) for i in z] # medians

    # # means
    # for xpos, ypos, yname in zip(x, y, ydescrip):
    #     plt.annotate(yname, (xpos, ypos), xytext=(0, 0), va='bottom',textcoords='offset points')
    # plt.annotate('Means (black)', (0.75, 50), xytext=(-10, 60), va='bottom',textcoords='offset points',arrowprops=dict(facecolor='black', shrink=0.05))

    # # medians (red)
    # for xpos, zpos, zname in zip(x, z, zdescrip):
    #     plt.annotate(zname, (xpos, zpos), xytext=(0, 0), va='bottom',textcoords='offset points',color='Red')
    # plt.annotate('Medians (red)', (0.75, 15), xytext=(-10, -60), va='bottom',textcoords='offset points',color='Red',arrowprops=dict(facecolor='Red', shrink=0.05))

    # # plotting
    # plt.xlabel('Collapse clades with average distance to leaves < X (e.g. 0.2)')
    # plt.ylabel('Mean size of clusters after collapse')
    # plt.show()


#-----------------------------------------------------------------#
#   @Example                                                       #
#-----------------------------------------------------------------#
#
#
# collapsed at distances,d, from d_from to d_to, with d_step increments
#
#
#
# @WARN: here is where you set the parameters
# @PARAMETERS    :   note that save_d determines the d-cutoff to save the cluster data at
#
d_from  = 0.0
d_to    = 0.3
d_step  = 0.002
e       = 0.05

d_vec = np.arange(d_from,d_to,d_step)

for d_cut in d_vec:

    
    #
    #
    save_d  = d_cut
    #
    #
    #
    print('###################')
    print('Distance Cut-off: '+str(save_d))
    print('###################')
    #
    #
    #
    #
    if not os.getcwd()=='/home/ab108/0VB/2kb/scripts/stamp_to_collapsed':
        print('You need to be in "./stamp_to_collapsed" ... attempting to get there now automatically ...')
        os.chdir('./stamp_to_collapsed/')
        if not os.getcwd()=='/home/ab108/0VB/2kb/scripts/stamp_to_collapsed':
            print('Failed to get to: "./stamp_to_collapsed", do so manually then re-run... ')
            print('\t Full Path: /home/ab108/0VB/2kb/scripts/stamp_to_collapsed ...')
    # try: 
    #     cut_tree = exploreCutoffs( '../../data/stamp_data/out/dreme_100bp_e'+str(e)+'/SWU_SSD/', 'out.tree', d_from = d_from, d_to = d_to, d_step = d_step,  save_d = save_d )
    # except: 

    #
    #   MAIN PROGRAM HERE >>>>>>
    #
    #                                       TREE FILE LOCATION                                TREE NAME
    #
    try:
        cut_tree = exploreCutoffs( '../../data/stamp_data/out/dreme_100bp_e'+str(e)+'/SWU_SSD/', 'out.tree', d_from = d_from, d_to = d_to, d_step = d_step,  save_d = save_d )

    except ete2.parser.newick.NewickError:
        print(' E.T.E. says: ete2.parser.newick.NewickError: Unexisting tree file or Malformed newick tree structure.')
        raise Exception('\tThe most common problem is that the tree is missing from: ../data/stamp_data/out/dreme_100bp_e0.05/SWU_SSD/out.tree')
    #
    #
    #
    #--------------------------#
    # Generate fasta per motif #
    #--------------------------#
    #
    #
    # Dictionary for dreme motifs names -> transfac pwms    ; 
    #
    print('Reading motifs from dreme output...')
    dataPath_dremeAllSp = '../../data/stamp_data/in/dreme_100bp_e'+str(e)+'/'          # clustered motif fastas, given d
    filename_dremeAllSp = 'dreme.fasta'
    file_in_dremeAllSp  = open(dataPath_dremeAllSp+filename_dremeAllSp,'r')

    print('Reading dreme.fasta for motif headers --> motif headers are then keyed to their FASTAs --> to later pump them into stamp')

    motif_to_pwm = {}

    while True:
        line = file_in_dremeAllSp.readline()
        if line == "":                      # break when finished
            break
        
        elif line.startswith('DE '):
            header = line.split(' ')[1].rstrip()
            motif_to_pwm[header] = []
            
            while True:
                line = file_in_dremeAllSp.readline()
                if "XX" in line:    # break when fini
                    break
                else:
                    motif_to_pwm[header].append(line)
    file_in_dremeAllSp.close()

    # Cluster file generator
    #
    dataPath_clusters_in    = '../../data/stamp_data/out/dreme_100bp_e'+str(e)+'/SWU_SSD/cluster_motifs_d'+str(save_d)+'/headers_for_fasta/'

    filename_clusters_list  = os.listdir(dataPath_clusters_in) # list the collapsed motif files to iterate over

    print(' Generating /fasta/ system... headers_for_fasta --> as input to pull out fastas ')

    dataPath_clusters_out   = '../../data/stamp_data/out/dreme_100bp_e'+str(e)+'/SWU_SSD/cluster_motifs_d'+str(save_d)+'/fasta/'

    # OVERWRITE /FASTA/ FOLDER
    if os.path.exists(dataPath_clusters_out):
        print('cluster fasta folder already exists, overwriting...')
        import shutil
        shutil.rmtree(dataPath_clusters_out)
        os.makedirs(dataPath_clusters_out)
    else:
        os.makedirs(dataPath_clusters_out)

    # Dictionary of clusters  -> motifs -> pwms 

    cluster_to_motif_to_pwm = {}

    for filename in filename_clusters_list:     # list of clusters

        print('\t'+filename)

        filename_fasta  = filename.replace('.txt','.fasta')
        clusterName     = filename.replace('.txt','')

        if not cluster_to_motif_to_pwm.has_key(clusterName):                    # had_key checks for duplicates, which we don't want!
            cluster_to_motif_to_pwm[clusterName] = {}                               
        else:
            raise Exception('WARNING: duplicate clustername? How is that possible..?')

        file_cluster_in = open(dataPath_clusters_in +filename,'r')              # file to read transfac motifs from dreme
        file_cluster_out= open(dataPath_clusters_out+filename_fasta,'w')        # file to write the clusters into

        while True:
            motif = file_cluster_in.readline().rstrip()
            if motif == "":
                break
            motif_transfac                              = motif_to_pwm[motif]                                    # get the transfac lines 
            cluster_to_motif_to_pwm[clusterName][motif] = motif_transfac            # store
            file_cluster_out.write('DE '+motif+'\n'+''.join(motif_transfac)+'XX\n') # write to the cluster file

        file_cluster_in.close() 
        file_cluster_out.close()


    #-----------------------------------------------------------------
    # STAMP
    #-----------------------------------------------------------------

    # path args
    dataPath_i_fasta        = '../../data/stamp_data/out/dreme_100bp_e'+str(e)+'/SWU_SSD/cluster_motifs_d'+str(save_d)+'/fasta/'
    dataPath_i_common       = '../../data/stamp_data/in/common/'
    dataPath_o              = '../../data/stamp_data/out/dreme_100bp_e'+str(e)+'/SWU_SSD/cluster_motifs_d'+str(save_d)+'/family_binding_matrices/'

    # OVERWRITE /FASTA/ FOLDER
    if os.path.exists(dataPath_o):
        print('family_binding_matrices folder already exists, overwriting...')
        import shutil
        shutil.rmtree(dataPath_o)
        os.makedirs(dataPath_o)
    else:
        os.makedirs(dataPath_o)

    # filenames
    filename_i_jasparDb     = 'jaspar.motifs'
    filename_i_scoreDist    = 'ScoreDists/JaspRand_SSD_SWU.scores'
    filename_i_fasta_list   = os.listdir(dataPath_i_fasta) # list the collapsed motif cluster files (transfec pmw fastas)

    # iterate stamp over collapsed motif clusters

    for filename_fasta in filename_i_fasta_list:

        print '----------------------------------------------------------'
        print filename_fasta # print motif names which stamp will process

        motif_out_path = dataPath_o+filename_fasta.replace('.fasta','') # path for outputs

        # COMMAND LINE: stamp -tf ./in/dreme.fasta -sd ./in/ScoreDists/JaspRand_SSD_SWU.scores -match ./in/jaspar.motifs -cc SSD -align SWU -out ./out/hehe

        stamp = [   'stamp',
                    '-tf',      str(dataPath_i_fasta +filename_fasta),      # verbosity, 1:5
                    '-sd',      str(dataPath_i_common+filename_i_scoreDist),# over-write directory with <name> and write in the results
                    '-match',   str(dataPath_i_common+filename_i_jasparDb), # database of motifs to match against
                    '-cc',      'SSD',                                      # distance calcualtion method
                    '-align',   'SWU',                                      # alignment method, e.g. smith-waterman ungapped alignment 
                    '-out',     motif_out_path                              # output file prefix to match motif names
                    ]
        all_stamp = [stamp]
        for stamp in all_stamp:
            a = datetime.datetime.now()
            subprocess.Popen(stamp,shell=False).communicate()               # call commands to the shell with messages
            #subprocess.Popen(stamp,shell=False)
            b = datetime.datetime.now()
            #print('\t'+str(b-a)) # print st

    print('COMPLETE!')
