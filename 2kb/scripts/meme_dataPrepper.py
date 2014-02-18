import subprocess

# Need to ensure data is ready for MEME

# TEST MODE:
# file_in_name    = '../data/meme_data/test/anopheles_gambiae_upstream_memeready_100.fasta'
# file_out_name   = '../data/meme_data/test/anopheles_gambiae_upstream_memeready_independent.fasta'

# REAL MODE:
file_in_name    = '../data/meme_data/test/anopheles_gambiae_upstream.fasta'
file_out_name   = '../data/meme_data/test/anopheles_gambiae_upstream_memeready_all.fasta'

file_in  = open(file_in_name)        # Read file in order to establish filters that will later format the data ready for MEME
file_out = open(file_out_name,'w')
headers = {}
print 'Establishing filters...'
while True:
    header  = file_in.readline()
    sequence= file_in.readline()
    if header == "":                    # break when finished
        break
    header_split = header.split('\t')
    geneId      = header_split[0].replace('>','')
    utrCoord    = header_split[1].replace('UtrCoord:','').split(':')
    utrLength   = header_split[2].replace('UtrLength:','')
    utrType     = header_split[3].replace('UtrType:','').rstrip()
    chromosome  = utrCoord[0]
    start       = utrCoord[1].split('-')[0]
    end         = utrCoord[1].split('-')[1]
    strand      = utrCoord[2]
    headers[geneId] = {'utrCoord':{'chromosome':chromosome,'start':start,'end':end,'strand':strand},'utrLength':utrLength,'utrType':utrType}
file_in.close()


# Pulley samples:       shared samples with genes in the same direction 
#
#   ---S--- + ---G-------->
#   ---S--- + ---G--->
#   ---S--- + ---G------------>     ...only one of these samples are kept, rest are thrown...
# 
#   dealt with by making a dictionary that ensures a 1 to 1 key-value relationshipof sample-to-geneId
#
#   during the next .readline() loop utrcoord-to-geneid value will be checked against the current geneId in current iteration
#
geneIds = headers.keys()
one_coord_per_geneId = {}
for geneId in geneIds:
    utrCoord = headers[geneId]['utrCoord']['chromosome']+':'+headers[geneId]['utrCoord']['start']+'-'+headers[geneId]['utrCoord']['end']+':'+headers[geneId]['utrCoord']['strand']
    one_coord_per_geneId[utrCoord] = geneId   # To ensure only one geneId = one utrCoord, I purposefully allow values to be overidden for a given key 

# Tug'o'war samples:    shared samples with genes in the opposite direction 
#
#   s1  <-------G--- + ---S---
#   s2                 ---S--- + ---S------->
#
#   dealt with by making a set of unique utrcoords stripped of their strand property, now upon encountering a match forward samples and reverse samples take their half of the tug-o-war sequence
#
coords          = [ ':'.join(i.split(':')[0:-1]) for i in one_coord_per_geneId.keys() ]     # remove strand property
coords_shared   = set([i for i in coords if coords.count(i)>1])

print 'Generating meme-ready data...'

file_in = open(file_in_name) # Open file for formatted parsing
collection = []
coords_shared_geneIds = []

print ('\tFiltering shared samples...')

while True:    
    header  = file_in.readline()
    sequence= file_in.readline().rstrip()
    if header == "": # break when finished
        break
# MEME filter:  =>  Is this sample seq too short for meme? 
    # YES:
    if len(sequence)<10:
        if sequence=='\n':
            collection.append(header)
        continue
    # NO:
    else:
        geneId  = header.split('\t')[0].replace('>','')
        utrCoord= header.split('\t')[1].replace('UtrCoord:','').split(':')
    # PULLEY filter:            =>  Is this GeneId unique? or permitted as a winning pulley sample?
        if one_coord_per_geneId[':'.join(utrCoord)] == geneId:
            utrCoord_noStrand   = ':'.join(utrCoord[0:1+1])
        # TUG'O'WAR filter:     =>  Is this sample in a tug'o'war? 
            if utrCoord_noStrand in coords_shared:
                utrCoord_onlyStrand = utrCoord[2]
                coords_shared_geneIds.append(header.split('\t')[0].replace('>',''))
                sequence = sequence[(len(sequence)/2):] # this works, since the rightermost end of the seq = gene, so we must sample from the end inwards
            # OKAY:             => then upon passing the filters, we can write the sample :)
            file_out.write(header+sequence+'\n')
# DUSTMASKER: low complexity region masker
print('\tMasking low complexity regions...')
subprocess.call(['dust',file_out_name,'10'],stdout=file_out) # TODO: what value do I use for 'cut-off', here given as 10

file_in.close()
file_out.close()


#================================================================================================
# PULLEY Samples:       exactly shared sample seqs, such that                   SOLVED
#================================================================================================
#
# s1    ---S--- + ---G------->
# s2    ---S--- + ---G--->
#
#   ...notice the two forward facing genes share the sample upstream, 
#   ...to prevent MEME from over-representing this we allocate it to just one of the samples
#
#
#================================================================================================
# OVERLAP sampels:      partially shared Sample Seqs                            UNSOLVED
#================================================================================================
#
#  Can be a problem to try and distirbute out the overlapping sample_seqs <==...
# 
# Samples: 
#  s1             ----------------------
#  s2       --------------------
#  s3                             -----------
#  s4         -----------
#  s5             ------------------
# 
#       ...well then how can we possibly distribute out this overlapping region between correponding samples? 
#          ...there seems only one option: remove the intersect! => but then our data will be so sparse!
#
# After rm intersect:
#
#  s1        |   ----------------------|
#  s2      --|------------------       |
#  s3        |                    -----|------
#  s4        |-----------              |
#  s5        |   ------------------    |
# 
#           ^ Not much data left!!!     ^
#
#
#================================================================================================
# TUG'O'WAR Samples:    samples that are shared by opposite facing genes        SOLVED
#================================================================================================
#
#
# Reap all coordinates that are fully shared by two adjacent genes which are in opposite orientation, a.k.a "Tug-Of-War" orientation
#
#   s1  <-------G--- + ---S---
#   s2                 ---S--- + ---S------->
#
#   ...notice the shared sample sequence between: 
#                                                 s1 (a reverse strand gene's upstream sample) and
#                                                 s2 (a forward strand gene's upstream sample)
#
#   ...well if this occurs, we give half to s1 and half to s2... 
#
#   s1                 ---S...
#   s2                 ...S---
#
#


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ================================================================
# OLD METHODS:
# ================================================================

# Original
# print 'reading file...'
# collection = []
# while True:
#     header  = file_in.readline()
#     sequence= file_in.readline()
#     if header == "":                    # break when finished
#         break
#     if len(sequence)<10:
#         if sequence=='\n':
#             collection.append(header)
#         continue
#     else:
#         file_out.write(header+sequence)