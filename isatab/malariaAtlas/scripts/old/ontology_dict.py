import itertools
import numpy as np 
import difflib
import os 

#dataPath_out= '../data/isatab/'
dataPath_in = '../data/ontology/'
suffixes_in = ['country_id','species']


# SPECIES:ONTOTERMS         :       
filename_in = dataPath_in+'species'
file_in     = open(filename_in)
species_to_ontoTerm = {}
while True:
    line = file_in.readline()
    if line == "":                  # break when finished
        break
    line_split  = line.split('|')
    try:
        sp_miro = ' '.join(line_split[2].split())
        id_miro = ' '.join(line_split[1].split())
    except IndexError:
        continue
    species_to_ontoTerm[sp_miro] = id_miro
file_in.close()

sp_miro = species_to_ontoTerm.keys()


#--------------------------------------------------------------------------------

# country:ONTOTERMS         :       
filename_in = dataPath_in+'country_id'
file_in     = open(filename_in)

country_to_ontoTerm = {}
while True:
    line = file_in.readline()
    if line == "":                  # break when finished
        break
    line_split  = line.split('|')
    try:
        country_id  = ' '.join(line_split[1].split())
        country     = ' '.join(line_split[3].split())
        gazTerm     = ' '.join(line_split[4].split())
        un          = ' '.join(line_split[2].split())
    except IndexError:
        continue
    country_to_ontoTerm[country_id] = {'country':country,'gazTerm':gazTerm,'un':un}
file_in.close()

#--------------------------------------------------------------------------------

# SP_RAW:SP_MIRO
dataPath_out= '../data/isatab/'
dataPath_in = '../data/raw/'
suffix_in   = '_DVS.txt'
continents  = ['Asia','Africa','Americas','Europe']

print('generating isatabs...')
data = {}
for continent in continents: #iterate through each continent's raw spreadsheet
    print('\t'+continent)
    if not os.path.exists(dataPath_out+continent):
        print('No directory for: '+continent+' ...creating...')
        os.makedirs(dataPath_out+continent)

# MANAGE RAW SPREADSHEETS

    # READ DATA
    print('\t\tReading raw data...')
    filename_in = dataPath_in+continent+suffix_in
    file_in     = open(filename_in,'rU') #'rU' seems to fix a bug where newlines arn't read properly
    lines = []
    while True:
        line = file_in.readline()
        if line == "":
            break
        lines.append([i.rstrip() for i in line.split('\t')])
    file_in.close()

    # PAD MISSING DATA      :   Europe and Americas are missing columns in latter rows => uneven array => we pad rows to ensure even array by len(row) = 17
    print('\t\tPadding missing data...')
    lines_even = []
    for row in lines:
        if len(row) < 17:
            row = row+['']*(17-len(row))
        lines_even.append(row)

    # VECTORIZE DATA
    print('\t\tVectorizing data...')
    lines_arr       = np.array(lines_even)
    lines_arr_T     = lines_arr.T
    data[continent]            = {}
    for col in lines_arr_T:
        data[continent][col[0]]= col[1:]
    sampleNames     = list(data[continent]['id'])
    sampleMethods1  = list(data[continent]['sample method1'])
    sampleMethods2  = list(data[continent]['sample method2'])
    sampleMethods3  = list(data[continent]['sample method3'])
    sampleMethods4  = list(data[continent]['sample method4'])
    speciesMethods1 = list(data[continent]['id method1'])
    speciesMethods2 = list(data[continent]['id method2'])
    speciesNames    = [species.replace('"','') for species in list(data[continent]['species'])] #replace '"' with ''   :   some species are encapsulated in double-quotes
    countryIds      = list(data[continent]['country id'])
    latitudes       = list(data[continent]['latitude'])
    longitudes      = list(data[continent]['longitude'])
    year_starts     = list(data[continent]['year start'])
    year_ends       = list(data[continent]['year end'])
    month_starts    = list(data[continent]['month start'])
    month_ends      = list(data[continent]['month end'])
    dates           = [year+'-'+month_starts[i].zfill(2)+'-00'+'/'+year_ends[i]+'-'+month_ends[i].zfill(2)+'-00' for i,year in enumerate(year_starts)]
    citations       = list(data[continent]['citation'])
    assi            = list(data[continent]['ASSI'])



species_raw_unique = list(set(list(itertools.chain(*[[species.replace('"','') for species in list(data[continent]['species'])] for continent in continents]))))
species_raw_to_miro = {}
for species_raw in species_raw_unique:
    if species_raw_to_miro.has_key(species_raw):
        print('wtf')
    species_miro_match                  = difflib.get_close_matches(species_raw,sp_miro)[0]
    species_raw_to_miro[species_raw]    = {'sp_miro':species_miro_match,'id_miro':species_to_ontoTerm[species_miro_match]}








