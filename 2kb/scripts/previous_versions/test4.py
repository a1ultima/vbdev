

import os 
home = os.path.expanduser('~')
os.chdir(home)


file_in = open('0VB/2kb/data/kamil/Insecta_OGs.txt')                    # Open file for reading


gene_to_go = {}


while True:
    line = file_in.readline()
    if line == "":                  # break when finished
        break
    line_split = line.split('\t')   # split line where '\t' occurs

    geneId  = line_split[0]
    go_name = line_split[1]
    go_id   = line_split[2]

    if gene_to_go.has_key(geneId):
        gene_to_go[geneId].append(go_name)
    else:
        gene_to_go[geneId] = [go_name]

file_in.close()