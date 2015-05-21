dataPath_in = '../data/isatab/'

continents = ['Africa','Americas','Europe','Asia']

print('Duplicate rows in ISA-Tabs shall be pruned...')
for continent in continents:
    print'\t'+continent
    filename_in_samples     = dataPath_in+continent+'/s_samples.txt'
    filename_in_collections = dataPath_in+continent+'/a_collections.txt'
    filename_in_species     = dataPath_in+continent+'/a_species.txt' 
    data = {'s_samples':{},'a_collections':{},'a_species':{}}    

# Samples
    # file_in_samples = open(filename_in_samples,'r') 
    # headers = file_in_samples.readline() # skip the headers
    # while True:
    #     line        = file_in_samples.readline()
    #     line_split  = line.rstrip().split('\t') # split line where '\t' occurs
    #     if line == "": # break when finished
    #         break
    #     sourceName  = line_split[0]
    #     sampleName  = line_split[1]
    #     description = line_split[2]
    #     material    = line_split[3]
    #     term        = line_split[4]
    #     accession   = line_split[5]
    #     primaryKey  = '\t'.join([description,material,term,accession])
    #     if data['s_samples'].has_key(primaryKey):
    #         data['s_samples'][primaryKey].append('\t'.join([sourceName,sampleName])+'\t')
    #     else: 
    #         data['s_samples'][primaryKey] = ['\t'.join([sourceName,sampleName])+'\t']
    # file_in_samples.close()

# Collections
    print '\t\t pruning: a_collections.txt...'
    # Reading
    file_in_collections = open(filename_in_collections,'r') 
    headers = file_in_collections.readline() # skip the headers
    while True:
        line        = file_in_collections.readline()
        line_split  = line.split('\t') # split line where '\t' occurs
        if line == "": # break when finished
            break
        sourceName  = line_split[0]
        assayName   = line_split[1]
        description = line_split[2]
        protocol1   = line_split[3]
        protocol2   = line_split[4]
        protocol3   = line_split[5]
        protocol4   = line_split[6]
        date        = line_split[7]
        country     = line_split[8]
        term        = line_split[9]
        accession   = line_split[10]
        latitude    = line_split[11]
        longitude   = line_split[12].rstrip()
        primaryKey  = '\t'.join([description,protocol1,protocol2,protocol3,protocol4,date,country,term,accession,latitude,longitude])
        if data['a_collections'].has_key(primaryKey):
            data['a_collections'][primaryKey].append('\t'.join([sourceName,assayName])+'\t')
        else: 
            data['a_collections'][primaryKey] = ['\t'.join([sourceName,assayName])+'\t']
    file_in_collections.close()
    # Writing
    filename_out_collections= dataPath_in+continent+'/a_collections_pruned.txt'
    file_out_collections    = open(filename_out_collections,'w')
    file_out_collections.write(headers)
    for primaryKey in data['a_collections'].keys():
        uniqueSample = data['a_collections'][primaryKey][0]
        file_out_collections.write(uniqueSample+primaryKey+'\n')
    file_out_collections.close()

# Species
    print '\t\t pruning: a_species.txt...'
    # Reading
    file_in_species = open(filename_in_species,'r') 
    headers = file_in_species.readline() # skip the headers
    while True:
        line        = file_in_species.readline()
        line_split  = line.split('\t') # split line where '\t' occurs
        if line == "": # break when finished
            break
        sourceName  = line_split[0]
        assayName   = line_split[1]
        description = line_split[2]
        protocol1   = line_split[3]
        protocol2   = line_split[4]
        species     = line_split[5]
        term        = line_split[6]
        accession   = line_split[7]
        assi        = line_split[8].rstrip()
        primaryKey  = '\t'.join([description,protocol1,protocol2,species,term,accession])
        if data['a_species'].has_key(primaryKey):
            data['a_species'][primaryKey].append('\t'.join([sourceName,assayName])+'\t')
        else: 
            data['a_species'][primaryKey] = ['\t'.join([sourceName,assayName])+'\t']
    file_in_species.close()
    # Writing
    filename_out_species= dataPath_in+continent+'/a_species_pruned.txt'
    file_out_species    = open(filename_out_species,'w')
    file_out_species.write(headers)
    for primaryKey in data['a_species'].keys():
        uniqueSample = data['a_species'][primaryKey][0]
        file_out_species.write(uniqueSample+primaryKey+'\n')
    file_out_species.close()