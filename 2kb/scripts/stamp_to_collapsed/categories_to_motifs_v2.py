
import os.path

home = os.path.expanduser('~')

os.chdir(home)

os.chdir('..')

merged_motifs_d0p002_clade3_in = 'maccallr/agcc/fimo/merged-motifs-d0.002-species3'

motif_files = os.listdir(merged_motifs_d0p002_clade3_in)

motif_dict = { 'shuffled':{}, 'observed':{} }



# Dict so that we can focus on just SHUFFLED or OBSERVED filenames

for motif_file in motif_files:
    
    if '.txt' in motif_file:
    
        if '.shuffled' in motif_file:
    
            motif_dict['shuffled'][motif_file] = ''
    
        else:
    
            motif_dict['observed'][motif_file] = ''



# Dict of filename -to- clade -to- motif name

motif_list = motif_dict['observed'].keys()



# motif -to- clade version, i.e. clades are collapsed 

m = {}

for fi_name in motif_list:

    fi = open(merged_motifs_d0p002_clade3_in+'/'+fi_name,'r')

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
#   - which file do i use as input? e.g./home/maccallr/agcc/fimo/merged-motifs-d0.002-clade3/ anopheles.meme.txt          dipteran.meme.txt            g_complex.meme.txt          mosquito.meme.txt
#   - 

fimo_anopheles_gambiae = 'maccallr/agcc/fimo/merged-motifs-d0.002-species3/all-concatenated/fimo_anopheles_gambiae/fimo.report-expr.txt'

file_fimo_anopheles_gambiae = open(fimo_anopheles_gambiae,'r')


#--------------------------------------------------------------------------------


# creates dict where clade points to index of the empty matrix

clade_list = list(set(m.values()))

# generate empty keys, where keys are clade...

#   (1) start with ONE_TUPLE clades...

print('initiating each clade counter key with 0s...')

one_tuple_clade_TO_count = {}

for one_tuple_clade in clade_list:

    one_tuple_clade_TO_count[one_tuple_clade] = 0

# read through file line-by-line, each line is a 3-tuple of motifs, each belonging to a given species...

print('reading exprmap_to_motif_trio file. counting no. rows with 1-tuple-species X...')

while True:

    line = file_fimo_anopheles_gambiae.readline()
    
    if line == '':
        break
    
    if 'motif tuple' in line:

        # list of space-delimited readlines 
        motifs = line.split(' ')

        # the three clades that are present in the row 
        clade1 = m[motifs[3]]; clade2 = m[motifs[5]]; clade3 = m[motifs[7]]

        # iteratively count +1 if the incoming tuple contrains the clade
        one_tuple_clade_TO_count[clade1] = one_tuple_clade_TO_count[clade1] + 1

        one_tuple_clade_TO_count[clade2] = one_tuple_clade_TO_count[clade2] + 1
        
        one_tuple_clade_TO_count[clade3] = one_tuple_clade_TO_count[clade3] + 1

file_fimo_anopheles_gambiae.close()

#   (2) TWO_TUPLE clades...


