from ete2 import Tree

# FUNCTIONS:

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
            avg_distance_to_tips = mean([root_distance[tip]-root_distance[node]
                                         for tip in node2tips[node]])

            if avg_distance_to_tips < min_dist:
                # do whatever, ete support node annotation, deletion, labeling, etc.

                # rename
                node.name += ' COLLAPSED avg_d:%g {%s}' %(avg_distance_to_tips,
                                                 ','.join([tip.name for tip in node2tips[node]]))
                # label
                node.add_features(collapsed=True)

                # set drawing attribute so they look collapsed when displayed with tree.show()
                node.img_style['draw_descendants'] = False

                # etc...

#------------------------------------------------------------------------------------------------------

# Load a tree structure from a newick file.
dataPath_in = '../../data/stamp_data/SWU_SSD/'
filename_in = 'out.tree'
t = Tree(dataPath_in+filename_in)


import numpy as np 

#d = 0.001
clusterSize_dist_arr= []
clusterSize_mean_arr= []
clusterSize_std_arr = []

distances = np.arange(0,1,0.05)

for count,d in enumerate(distances):     # explore distance threshold params
    if count%10 == 0:
        print(count)
# apply collapse
    t = Tree(dataPath_in+filename_in)           # load tree
    collapse(t, d)                              # collapse /u threshold d
    for n in t.search_nodes(collapsed=True):    # re-annotate via collapses
        for ch in n.get_children():
            ch.detach()
# node distances
    dist_dict = cache_distances(t)
# unite nested list of collapsed motifs 
    clusters_multiple   = [i.name.split('{')[1][:-1].split(',') for i in dist_dict.keys() if "COLLAPSED" in i.name] # list of clusters, where each cluster is a list of motifs)
    clusters_single     = [[i.name] for i in dist_dict.keys() if not 'NoName' in i.name]
    clusters_all        = clusters_single+clusters_multiple
# ranked distribution of cluster sizes, largest-to-smallest
    clusterSize_dist= np.array(sorted([float(len(i)) for i in clusters_all],reverse=True))
    clusterSize_mean= np.mean(clusterSize_dist) 
    clusterSize_std = np.std(clusterSize_dist)
# store
    clusterSize_dist_arr.append(clusterSize_dist) 
    clusterSize_mean_arr.append(clusterSize_mean)
    clusterSize_std_arr.append(clusterSize_std)

import matplotlib.pyplot as plt
x = distances
y = clusterSize_mean_arr
yerr = [0]*len(distances)
xerr = clusterSize_std_arr
descrip = [str(int(round(i))) for i in y]
plt.errorbar(x, y, xerr, yerr, capsize=0, ls='none', color='black', elinewidth=1)
for xpos, ypos, name in zip(x, y, descrip):
    plt.annotate(name, (xpos, ypos), xytext=(0, 0), va='bottom',textcoords='offset points')
plt.show()

print(clusterSize_mean_arr)

clusterStats = [clusterSize_dist_arr,clusterSize_mean_arr,clusterSize_std_arr]

# import pickle
# import os.path
# pickle_name = 'clusterStats'
# pickle_obj  = clusterStats
# if os.path.isfile(pickle_name+'.p'):
#     pickle_obj = pickle.load(open(pickle_name+'.p','rb'))
# else:
#     pickle.dump(pickle_obj,open(pickle_name+'.p','wb'))