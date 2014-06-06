import subprocess
import os

# Need to ensure data is ready for MEME

fileName_in = '../data/meme_data/test/anopheles_gambiae_upstream.fasta'
fileName_out= '../data/meme_data/test/anopheles_gambiae_upstream_memeready_all.fasta'
file_in     = open(fileName_in)                   # Open file for reading
file_out    = open(fileName_out,'w')

headers = {}

print 'reading file...'
while True:
    header  = file_in.readline()
    sequence= file_in.readline()
    if header == "":                    # break when finished
        break

    header_split = header.split('\t')

    geneId      = header_split[0]
    utrLength   = int(header_split[2].replace('UtrLength:',''))
    utrType     = header_split[3].replace('UtrType:','')

    utrCoord    = header_split[1].replace('UtrCoord:','').split(':')

    chromosome  = utrCoord[0]
    start       = int(utrCoord[1].split('-')[0])
    end         = int(utrCoord[1].split('-')[1])
    strand      = utrCoord[2]

    headers[geneId] = {'utrCoord':{'chromosome':chromosome,'start':start,'end':end,'strand':strand},'utrLength':utrLength,'utrType':utrType}

# Run Dustmasker for the low complexity regions

# os.chdir('../data/meme_data/test/')
# fileName_dust =fileName_out.split('/')[-1]
# subprocess.call(['dust',fileName_dust,'10'],stdout=file_out) # TODO: what value do I use for 'cut-off', here given as 10

file_in.close()
file_out.close()