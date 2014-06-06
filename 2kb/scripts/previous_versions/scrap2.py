
# moo = 'aasdkhgfgdfdfhgkudfhjkdfhjkdfhgjkh'

# moo_listed = list(moo)

# regions = [(0,2),(1,2),(3,4)]
# returner = [moo_listed[j[0]:j[1]+1] for j in regions]


# #[moo_listed[j[0]:j[1]+1] = 'X'*(j[1]-j[0]) for j in regions]

# print('moo')

# for region in regions:
#     print(moo_listed)
#     print(moo_listed[region[0]:region[1]+1])
#     if masking = 'hard':
#         moo_listed[region[0]:region[1]+1] = 'X'*(len(moo_listed[region[0]:region[1]+1]))
#     if masking = 'soft':
#         moo_listed[region[0]:region[1]+1] = [letter.lower for letter in moo_listed[region[0]:region[1]+1])]


# masking = 'hard'

# # Iterate masking:
# for region in repeats_regions:
#     #print(region)
#     # Are repeat.Location's start or end out of bounds? --> restrict them to sample_locations limits if so
#     if region[0]<0:
#         region[0]=0
#     if region[1]>len(sample_seq):     # TODO: is len() the actual matcher? 
#         region[1]=len(sample_seq)
#     #print(region)
#     # Replace all repeating sequence by appropriate masks
#     if masking=='hard':
#         sample_seq[region[0]:region[1]+1]='N'*(len(sample_seq[region[0]:region[1]+1])) # PyCogent is not happy with X
#     elif masking=='soft':
#         sample_seq[region[0]:region[1]+1]=[letter.lower() for letter in sample_seq[region[0]:region[1]+1]]
# moo = ''.join(sample_seq)


"""

meme anopheles_gambiae_upstream_memeready_100.fasta -dna -mod anr -revcomp -minw 10 -maxw 12 -maxsize 500000

"""

geneIds = headers.keys()

#coords = [':'.join([headers[i]['utrCoord']['chromosome'],headers[i]['utrCoord']['start'],headers[i]['utrCoord']['end']) for i in geneIds]


coords  = [  headers[i]['utrCoord']['chromosome']+':'+headers[i]['utrCoord']['start']+'-'+headers[i]['utrCoord']['end'] for i in geneIds ]

coords_shared = set([x for x in coords if coords.count(x) > 1])

