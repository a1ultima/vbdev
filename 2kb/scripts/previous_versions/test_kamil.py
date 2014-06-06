import os 
home = os.path.expanduser('~')
os.chdir(home)

def get_gambiae_genes(file_in,file_out):
    geneIds_gambiae = []
    while True:
        line = file_in.readline()
        if line == "":                  # break when finished
            break
        line_split = line.split('\t')   # split line where '\t' occurs
        geneId = line_split[3]
        if 'AGAP' in geneId:
            geneIds_gambiae.append(geneId)
            file_out.write(geneId+'\n')
    file_in.close()
    file_out.close()
    return(geneIds_gambiae)

def gene2Go_dictGen(file_in):
    geneId_to_go = {}
    line = file_in.readline() # HEADERS
    while True:
        line = file_in.readline()
        if line == "":                    # break when finished
            break
        line_split = line.split('\t')     # split line where '\t' occurs
        geneId  = line_split[0]
        go_name = line_split[1]
        go_id   = line_split[2]
        if geneId_to_go.has_key(geneId):
            geneId_to_go[geneId].append(go_id)
        else:
            geneId_to_go[geneId] = [go_id]
    file_in.close()
    return geneId_to_go

def Og2Gene_dictGen(file_in):
    og_to_geneId = {}
    line = file_in.readline() # Skip HEADERS
    while True:
        line = file_in.readline()
        if line == "":                  # break when finished
            break
        line_split = line.split('\t')   # split line where '\t' occurs
        orthologous_group   = line_split[0]
        geneId              = line_split[3]
        if 'AGAP' in geneId:
            if og_to_geneId.has_key(orthologous_group):
                og_to_geneId[orthologous_group].append(geneId)
            else:
                og_to_geneId[orthologous_group] = [geneId]
    file_in.close()
    return(og_to_geneId)

def get_geneIds_fromFasta(file_in,onlyRepeats):
    """
    Notes:

    Args:
        file_in = open('0VB/2kb/data/meme_data/in/anopheles_gambiae_upstream_memeready_all_dusted.fasta')                   # Open file for reading
    """
    file_in = open(file_in)                   # Open file for reading
    tg_query    = re.compile(r'(TG)\1{3,}')
    tg_samples  = []
    while True:
        header  = file_in.readline()
        seq     = file_in.readline()
        if header == "":                    # break when finished
            break
        if onlyRepeats:
            if tg_query.search(seq):
                tg_samples.append(header.split('\t')[0].replace('>',''))
        else:
            if not tg_quesry.search(seq):
                    
    file_in.close()
    return tg_samples

#-------------------------------------------------------------------------------------------
# RUN: 
#-------------------------------------------------------------------------------------------


# # Kamil

# #Get gambiae gene Ids:
# file_in = open('0VB/2kb/data/kamil/Insecta_OGs.txt')                    # Open file for reading
# file_out= open('0VB/2kb/data/kamil/gambiae_geneIds.txt','w')                    # Open file for reading
# geneIds_gambiae = get_gambiae_genes(file_in,file_out)

# #Generate dictionary for geneIds to GO terms annotated by VB Biomart
# file_in = open('/home/ab108/0VB/2kb/data/kamil/mart_export.txt')                     # Open file for reading
# geneId_to_go = gene2Go_dictGen(file_in)

# #Generate dictionary for gene family to geneIds
# file_in = open('0VB/2kb/data/kamil/Insecta_OGs.txt')                    # Open file for reading
# og_to_geneId = Og2Gene_dictGen(file_in)


# TGTG repeats 

# get gene_to_go
file_in = open('/home/ab108/0VB/2kb/data/tg_repeats/tg3_biomart_go.txt')                     # Open file for reading
geneId_to_go_tgRepeat = gene2Go_dictGen(file_in)

import itertools

goIds_all = [i for i in list(itertools.chain(*geneId_to_go_tgRepeat.values())) if i]

file_out = open('/home/ab108/0VB/2kb/data/tg_repeats/goIds_lumped.txt','w')                     # Open file for reading

for i in goIds_all:
    file_out.write(i+'\n')
file_out.close()
