from ete2 import Tree

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

# Example                
#t = Tree("((A,(B:0.1,C:0.1)i1:0.1)i2:0.5,((D:0.1,(E:0.1,F:0.1)i3:0.1)i4:0.5,G)i5:0.3);", format=1)
#                  /-1.0, A
#           /0.5, i2
#          |      |       /-0.1, B
#          |       \0.1, i1
#          |              \-0.1, C
#-1.0, NoName
#          |              /-0.1, D
#          |       /0.5, i4
#          |      |      |       /-0.1, E
#           \0.3, i5      \0.1, i3
#                 |              \-0.1, F
#                 |
#                  \-1.0, G

print '.'

# Load a tree structure from a newick file.
dataPath_in = '../../data/stamp_data/SWU_SSD/'
filename_in = 'out.tree'
t = Tree(dataPath_in+filename_in)

print t.write()

#print t.get_ascii(attributes=["dist", 'name'])

# Now we collapse
collapse(t, 0.01)

# print ascii diagram
#print t.get_ascii(attributes=['name'])
#                                       /-A
#      /i2 COLLAPSED avg_d:0.466667 {C,A,B}
#     |                                  |                            /-B
#     |                                   \i1 COLLAPSED avg_d:0.1 {C,B}
#     |                                                               \-C
#-NoName
#     |                                      /-D
#     |   /i4 COLLAPSED avg_d:0.166667 {F,D,E}
#     |  |                                  |                            /-E
#      \i5                                   \i3 COLLAPSED avg_d:0.1 {F,E}
#        |                                                               \-F
#        |
#         \-G


# interactive tree visualization will hide collapsed nodes -- zoom, select, etc
#t.show()

# collapsed nodes are labeled, so you locate them and prune them
for n in t.search_nodes(collapsed=True):
    for ch in n.get_children():
        ch.detach()
#print t
#   /-i2 COLLAPSED avg_d:0.466667 {C,A,B}
#--|
#  |   /-i4 COLLAPSED avg_d:0.166667 {F,D,E}
#   \-|
#      \-G

print str(t.write())
# and the the pruned newick
#(i2 COLLAPSED avg_d_0.466667 {C_A_B}:0.5,(i4 COLLAPSED avg_d_0.166667 {F_D_E}:0.5,G:1)1:0.3);


# dist_dict = cache_distances(t)
# #--------------------------------------------\
# import pickle
# import os.path
# pickle_name = 'dist_dict'
# pickle_obj  = dist_dict 

# if os.path.isfile(pickle_name+'.p'):
#     pickle_obj = pickle.load(open(pickle_name+'.p','rb'))
# else:
#     pickle.dump(pickle_obj,open(pickle_name+'.p','wb'))
# #--------------------------------------------/
# # get names of the collapsed nodes 
# clusters_multiple   = [i.name.split('{')[1][:-1].split(',') for i in dist_dict.keys() if "COLLAPSED" in i.name] # list of clusters, where each cluster is a list of motifs)
# clusters_single     = [[i.name] for i in dist_dict.keys() if not 'NoName' in i.name]
# clusters_all        = clusters_single+clusters_multiple


