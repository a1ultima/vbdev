import numpy as np 
import pickle
import os 

# Ensure we are in ./scripts/
if not os.getcwd().endswith('scripts'):
    os.chdir('../')
    print os.getcwd()

# Load the pickle file containing the motif stats for every parameterisation event generated in the latest execution of "motifStatistics...py"
e_cut               = str(raw_input('what is the e-value of the data you want to work with? e.g. 0.05'))
print('Reading pickled cluster statistics...')
e_cut_dir           = '../data/stamp_data/out/dreme_100bp_e'+e_cut+'/SWU_SSD/'
stats_file          = 'cluster_to_stats.p'
cluster_to_stats    = pickle.load(open(e_cut_dir+stats_file,'rb'))

# Dictionary for pointing d-cutoffs to line-of-best-fit statistics y = mx + c 
d_to_line = {}

# Write the stats into cluster_motifs_d<...>
for event in cluster_to_stats.keys():
    print('Tree Parameterisation: '+event)

    # Pair cluster's number of unique species vs. cluster's average entropy into a vector, nSp_vs_H, delimited with tabs
    print('\tpairing species diversity (conservedness) to mean entropy of each cluster...')

    clusters    = cluster_to_stats[event].keys()

    nSp         = [cluster_to_stats[event][i]['species']['unique']['n_unique'] for i in cluster_to_stats[event].keys()]
    
    H           = [cluster_to_stats[event][i]['cluster']['H'] for i in cluster_to_stats[event].keys()]
    
    nSp_vs_H    = zip(clusters,nSp,H)
    nSp_vs_H    = '\n'.join([str(i)+'\t'+str(j)+'\t'+str(k) for i,j,k in nSp_vs_H])

    # Name of file to write data to
    print('\tWriting data...')
    d           = event.replace('e'+str(e_cut)+'_d', '')
    cluster_dir = e_cut_dir+'cluster_motifs_d'+d+'/statistics.txt'

    # Line of best fit
    vs = np.array(zip(nSp,H))
    line = np.polyfit(vs[:,0], vs[:,1], 1)
    d_to_line[d] = {'L':line,'V':vs}

    # Write data
    file_out = open(cluster_dir,'wb')
    file_out.write(nSp_vs_H)
    file_out.close()

print('COMPLETE!')







