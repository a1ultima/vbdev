
import os.path

home = os.path.expanduser('~')
os.chdir(home)

merged_motifs_d0p002_species3_in = '/home/maccallr/agcc/fimo/merged-motifs-d0.002-species3'

motif_files = os.listdir(merged_motifs_d0p002_species3_in)
motif_dict = {'shuffled':{},'observed':{}}

# Dict so that we can focus on just SHUFFLED or OBSERVED filenames
for motif_file in motif_files:
    if '.txt' in motif_file:
        if '.shuffled' in motif_file:
            motif_dict['shuffled'][motif_file] = ''
        else:
            motif_dict['observed'][motif_file] = ''

# Dict of filename -to- clade -to- motif name

motif_list = motif_dict['observed'].keys()

# filename -to- motif -to- clade version
# m = {}
# for fi_name in motif_list:
#     fi = open(merged_motifs_d0p002_species3_in+'/'+fi_name,'r')
#     clade = fi_name.replace('.meme.txt','')  
#     m[fi_name] = {}
#     #m[fi_name][clade] = []
#     while True: 
#         line = fi.readline()
#         #print line
#         if line=='':
#             break
#         if 'MOTIF' in line:
#             m[fi_name][line.rstrip().replace('MOTIF ','')] = clade
#     fi.close()

# motif -to- clade version, i.e. clades are collapsed 
m = {}
for fi_name in motif_list:
    fi = open(merged_motifs_d0p002_species3_in+'/'+fi_name,'r')
    clade = fi_name.replace('.meme.txt','')  
    while True: 
        line = fi.readline()
        if line=='':
            break
        if 'MOTIF' in line:
            m[line.rstrip().replace('MOTIF ','')] = clade
    fi.close()

# fimo.report-expr.txt 

# Questions for Bob
#   - which file do i use as input? e.g./home/maccallr/agcc/fimo/merged-motifs-d0.002-species3/ anopheles.meme.txt          dipteran.meme.txt            g_complex.meme.txt          mosquito.meme.txt
#   - 

fimo_anopheles_gambiae = '/home/maccallr/agcc/fimo/merged-motifs-d0.002-species3/all-concatenated/fimo_anopheles_gambiae/fimo.report-expr.txt'

file_fimo_anopheles_gambiae = open(fimo_anopheles_gambiae,'r')


#--------------------------------------------------------------------------------


# creates dict where species points to index of the empty matrix

print( 'generating species to list index pointer... ' )

os.path.expanduser('home/ab108/')

if not os.getcwd() == '/home/ab108/0VB/2kb/scripts':
    os.chdir('./0VB/2kb/scripts/')

import speciesManage

species_list = speciesManage.generate_list('species_list.txt') 

species_to_index = dict() # species tuples -to- species distribution

count = 0

for species in species_list:

    species_to_index[species] = count

    count = count + 1



#-------------

import numpy as np 



species_matrix_2tuples = np.array([[]]*len(species_list))


while True:

    line = file_fimo_anopheles_gambiae.readline()
    
    if line == '':
        break
    
    if 'motif tuple' in line:
        motifs = line.split(' ')
        
        # motifs
        motif1   = motifs[3] # pull out the 
        motif2   = motifs[5]
        motif3   = motifs[7]

        # species 
        species1 = m[motif1]
        species2 = m[motif2]
        species3 = m[motif3]

        if 

        # count the species present
        species_list = (species1,species2,species3)
        
        species_dict[species_list] = dict(zip(species_list,map(species_list.count,species_list)))

file_fimo_anopheles_gambiae.close()


