
#coords_exact          = [  headers[i]['utrCoord']['chromosome']+':'+headers[i]['utrCoord']['start']+'-'+headers[i]['utrCoord']['end']+':'+headers[i]['utrCoord']['strand'] for i in geneIds ]

#coords_shared   = set([i for i in coords if coords.count(i)>1])

one_coord_per_geneId = {}

for geneId in geneIds:

    utrCoord = headers[geneId]['utrCoord']['chromosome']+':'+headers[geneId]['utrCoord']['start']+'-'+headers[geneId]['utrCoord']['end']+':'+headers[geneId]['utrCoord']['strand']

    # if coords_exact_to_geneIds.has_key(utrCoord):
    #     coords_exact_to_geneIds[utrCoord].append(geneId)
    # else:
    #     coords_exact_to_geneIds[utrCoord] = [geneId]
    
    # I purposefully allow values to be overidden for a given key to ensure that only one geneId = one utrCoord
    one_coord_per_geneId[utrCoord] = [geneId]

#utrCoords = one_coord_per_geneId.keys()
