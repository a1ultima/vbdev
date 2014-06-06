
from ete2 import Tree

# Loads a tree structure from a newick string. The returned variable ’t’ is the root node for the tree.
# t = Tree("(A:1,(B:1,(E:1,D:1):0.5):0.5);" )

# Load a tree structure from a newick file.
dataPath_in = '../../data/stamp_data/SWU_SSD/'
filename_in = 'out.tree'

t = Tree(dataPath_in+filename_in)

# We can also write into a file
# t.write(format=1, outfile="new_tree.nw")

# print
# print t.write()

def collapsed_leaf(node):
    if len(node2labels[node]) <= 0.01:
       return True
    else:
       return False

#t = Tree("((((a,a,a)a,a)aa, (b,b)b)ab, (c, (d,d)d)cd);", format=1)
#print t

# We create a cache with every node content
#node2labels = t.get_cached_content(store_attr="name")
#print t.write(is_leaf_fn=collapsed_leaf)

# We can even load the collapsed version as a new tree
#t2 = Tree( t.write(is_leaf_fn=collapsed_leaf) )
#print t2

#t.write(format=1, outfile="mooTree.txt")