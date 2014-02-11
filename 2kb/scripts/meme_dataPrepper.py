

# Need to ensure data is ready for MEME

# HEADER MODE:
# file_in_name    = '../data/meme_data/test/anopheles_gambiae_upstream_memeready_100.fasta'
# file_out_name   = '../data/meme_data/test/anopheles_gambiae_upstream_memeready_independent.fasta'

# REAL MODE:
file_in_name    = '../data/meme_data/test/anopheles_gambiae_upstream.fasta'
file_out_name   = '../data/meme_data/test/anopheles_gambiae_upstream_memeready_all.fasta'

file_in     = open(file_in_name)                   # Open file for reading
file_out    = open(file_out_name,'w')

headers = {}
print 'reading file...'
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
file_in.close()

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
geneIds         = headers.keys()
coords          = [  headers[i]['utrCoord']['chromosome']+':'+headers[i]['utrCoord']['start']+'-'+headers[i]['utrCoord']['end'] for i in geneIds ]
coords_shared   = set([i for i in coords if coords.count(i)>1])

file_in         = open(file_in_name)                   # Open file for reading again

print 'reading file...'
collection = []
count = 0
while True:
    count = count + 1
    header  = file_in.readline()
    sequence= file_in.readline().rstrip()
    if header == "":                    # break when finished
        break
    if len(sequence)<10:
        if sequence=='\n':
            collection.append(header)
        continue
    else:
        utrCoord            = header.split('\t')[1].replace('UtrCoord:','').split(':')
        utrCoord_noStrand   = ':'.join(utrCoord[0:1+1])

        if utrCoord_noStrand in coords_shared:
            utrCoord_onlyStrand = utrCoord[2]

            # ORIENTATION OPTION 1 :                                ----v   which of these options is correct? / Hint: think about the revcomp as seqs exit the db
            sequence = sequence[(len(sequence)/2):] # this works, since the rightermost end of the seq = gene, so we must sample from the end in
            # sequence = sequence[0:(len(sequence)/2)] # this fails

            # ORIENTATION OPTION 2 :                                ----^
            # if utrCoord_onlyStrand=='1':
            #     print('plus')
            #     sequence = sequence[(len(sequence)/2):]
            # elif utrCoord_onlyStrand=='-1':
            #     print('minus')
            #     sequence = sequence[0:(len(sequence)/2)+1]

        file_out.write(header+sequence+'\n')

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

file_in.close()
file_out.close()