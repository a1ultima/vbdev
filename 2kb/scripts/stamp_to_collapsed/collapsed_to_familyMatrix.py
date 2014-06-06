import datetime
import os 
import subprocess

save_d = 0.5

# path args
dataPath_i_fasta        = '../../data/stamp_data/out/dreme_100bp_e0.5/SWU_SSD/cluster_motifs_d'+str(save_d)+'/fasta/'
dataPath_i_common       = '../../data/stamp_data/in/common/'
dataPath_o              = '../../data/stamp_data/out/dreme_100bp_e0.5/SWU_SSD/cluster_motifs_d'+str(save_d)+'/family_binding_matrices/'
if not os.path.exists(dataPath_o):
    os.makedirs(dataPath_o)
# filenames
filename_i_jasparDb     = 'jaspar.motifs'
filename_i_scoreDist    = 'ScoreDists/JaspRand_SSD_SWU.scores'
filename_i_fasta_list   = os.listdir(dataPath_i_fasta) # list the collapsed motif cluster files (transfec pmw fastas)

# iterate stamp over collapsed motif clusters
for filename_fasta in filename_i_fasta_list:
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
        b = datetime.datetime.now()
        print('\t'+str(b-a)) # print duration

