import speciesManage

speciesFile         = './species_list.txt'
dataPath_in         = '../data/meme_data/out/dreme_100bp/sampled_all_hpc/'
suffix_in           = '_100bp/dreme.txt'
species_to_motifs   = {}
species_list        = speciesManage.generate_list(speciesFile)

# READ DREME => species:motifs
for species in species_list:
    print(species)
    species_to_motifs[species]  = []
    filename_in                 = dataPath_in+species+suffix_in
    file_in                     = open(filename_in,'r')
    count                       = 0
    while True:
        line = file_in.readline()
        if line == "":
            break
        if 'MOTIF' in line:
            count += 1
            motif = line.split(' ')[1].rstrip()   # split line where '\t' occurs
            species_to_motifs[species].append(motif)
    file_in.close()

