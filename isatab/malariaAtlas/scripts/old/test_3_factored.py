import itertools
import numpy as np
import os 
import re 
import difflib

# Pipeline:
#    1. Generate dictionaries for: species, country, citations  --> species_to_miro, country_to_gaz, citations_to_pubmedId
#    2. Re-interpret the above dictionaries to allow conversion from raw_to_<aboveDictionary>  --> e.g. raw_to_species_to_miro
#    3. Read & Vectorize raw data columns
#    4. Parse vectorized columns partitioned into: s_samples, a_collections, a_species

print('')
print('#######################')
print('#  ISATAB GENERATOR:  #')
print('#######################')
print('')

#-------------------------------------------------------------------------------------------- 

#===============#
# ONTOLOGY GEN  #
#===============#
print('Reading ontology dictionaries...')
dataPath_in = '../data/ontology/'
suffixes_in = ['country_id','species']

# SPECIES:ONTOTERMS         :  
print('\tSPECIESNAMES_to_MIROTERMS')     
filename_in = dataPath_in+'species'
file_in     = open(filename_in)
species_to_ontoTerm = {}
while True:
    line = file_in.readline()
    if line == "":                  # break when finished
        break
    line_split  = line.split('|')
    try:
        sp   = ' '.join(line_split[2].split())
        onto = ' '.join(line_split[1].split())
    except IndexError:
        continue
    species_to_ontoTerm[sp] = onto
file_in.close()
sp_miro = species_to_ontoTerm.keys()

# COUNTRY_ID:ONTOTERMS         :
print('\tCOUNTRIES_to_GAZTERMS')            
filename_in = dataPath_in+'country_id'
file_in     = open(filename_in)
country_to_ontoTerm = {}
while True:
    line        = file_in.readline()
    line_split  = line.split('|')
    if line == "":                  # break when finished
        break
    try:
        country_id  = ' '.join(line_split[1].split())
        country     = ' '.join(line_split[3].split())
        gazTerm     = ' '.join(line_split[4].split())
        un          = ' '.join(line_split[2].split())
    except IndexError:
        continue
    country_to_ontoTerm[country_id] = {'country':country,'gazTerm':gazTerm,'un':un}
file_in.close()

#-------------------------------------------------------------------------------------------- 

#===============#
#  ISATAB GEN   #
#===============#

dataPath_out= '../data/isatab/'
dataPath_in = '../data/raw/'
suffix_in   = '_DVS.txt'
continents  = ['Asia','Africa','Americas','Europe']

print('Formatting raw data...')
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
    data[continent] = {}
    for col in lines_arr_T:
        data[continent][col[0]]= col[1:]
    #speciesNames = [species.replace('"','') for species in list(data[continent]['species'])] #replace '"' with ''   :   some species are encapsulated in double-quotes

# SPECIES_RAW:(SPECIES_MIRO,ID_MIRO)
print('Generating raw_species:miro_species dictionary...')
species_raw_unique  = list(set(list(itertools.chain(*[[species.replace('"','') for species in list(data[continent]['species'])] for continent in continents]))))
species_raw_to_miro = {}
for species_raw in species_raw_unique:
    species_miro_match                  = difflib.get_close_matches(species_raw,sp_miro)[0]
    species_raw_to_miro[species_raw]    = {'sp_miro':species_miro_match,'id_miro':species_to_ontoTerm[species_miro_match]}

# WRITE ISATABS      :       s_samples, a_collections, a_species
print('writing ISA-TABs to file...')
for continent in continents:
    print('\t'+continent)
    # READ VECTORIZED DATA      :       Columns of data that will be divided into: s_samples, a_collections, a_species
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

    # I_INVESTIGATIONS HELPER   :   finds all unique species ID & sampling methods, it accumulate across continents
    try:
        print('\t\tAccumulating methods: i_investigations...')
        unique_species_methods = list(set(list(data[continent]['id method1'])+list(data[continent]['id method1'])+unique_species_methods))
        unique_sample_methods  = list(set(list(data[continent]['sample method1'])+list(data[continent]['sample method2'])+list(data[continent]['sample method3'])+list(data[continent]['sample method4'])+unique_sample_methods))
    except NameError:
        print('\t\tAccumulating methods: i_investigations...')
        unique_species_methods = list(set(list(data[continent]['id method1'])+list(data[continent]['id method1'])))
        unique_sample_methods  = list(set(list(data[continent]['sample method1'])+list(data[continent]['sample method2'])+list(data[continent]['sample method3'])+list(data[continent]['sample method4'])))

#===============#
#   S_SAMPLES   #
#===============#
    print('\t\tGenerating Sheets...')
    print('\t\t\tS_SAMPLES')
    suffix_out  = '/s_samples.txt'
    filename_out= dataPath_out+continent+suffix_out
    # ISA-TAB COLUMNS: ------------------------------------------------------------------------------------------------------------------------------
    #                       |                       |                   |               |               |                   |                       |
    headers     =           'Source Name\t          Sample Name\t       Description\t   Material Type\t Term Source REF\t   Term Accession Number\n'
    values      = ''.join([ 'malariaAtlas'+'\t'+    sampleName+'\t'+    ''+'\t'+        ''+'\t'         ''+'\t'+            ''+'\n'                 for sampleName in sampleNames]) 
    #                       |                       |                   |               |               |                   |                       |
    #------------------------------------------------------------------------------------------------------------------------------------------------
    headers     = re.sub(r'\t ','\t',re.sub(r'( )\1{2,}', ' ', headers))
    parse       = headers+values
    file_out    = open(filename_out,'w')
    file_out.write(parse)
    file_out.close()

#===================#
#   A_COLLECTIONS   #
#===================#
    print('\t\t\tA_COLLECTIONS')
    suffix_out  = '/a_collections.txt'
    filename_out= dataPath_out+continent+suffix_out
    # ISA-TAB COLUMNS: --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #                       |                   |                               |               |                       |           |                       |           |                       |           |                       |           |               |                                                   |                   |                                                                   |                                                           |                                                           |
    headers     =           'Sample Name\t      Assay Name\t                    Description\t   Protocol REF\t          Performer\t Protocol REF\t          Performer\t Protocol REF\t          Performer\t Protocol REF\t          Performer\t Date\t          Country\t                                           Term Source REF\t   Term Accession Number\t                                             Characteristics [Collection site latitude (VBcv:0000817)\t  Characteristics [Collection site longitude (VBcv:0000816)\n'
    values      = ''.join([ sampleName+'\t'+    sampleName+' collection'+'\t'+  ''+'\t'+        sampleMethods1[i]+'\t'+ ''+'\t'+    sampleMethods2[i]+'\t'+ ''+'\t'+    sampleMethods3[i]+'\t'+ ''+'\t'+    sampleMethods4[i]+'\t'+ ''+'\t'+    dates[i]+'\t'+  country_to_ontoTerm[countryIds[i]]['country']+'\t'+ 'GAZ'+'\t'+         country_to_ontoTerm[countryIds[i]]['gazTerm'].split(':')[1]+'\t'+   latitudes[i]+'\t'+                                          longitudes[i]+'\n'                                           for i,sampleName in enumerate(sampleNames)]) # tab-delim text to parse
    #                       |                   |                               |               |                       |           |                       |           |                       |           |                       |           |               |                                                   |                   |                                                                   |                                                           |                                                           |
    #--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #                                                                                                     ^--------------------^-----------------------^----------------------^ --> needs to be dictified 
    headers     = re.sub(r'\t ','\t',re.sub(r'( )\1{2,}', ' ', headers))
    parse       = headers+values
    file_out    = open(filename_out,'w')
    file_out.write(parse)
    file_out.close()

#===============#
#   A_SPECIES   #
#===============#
    print('\t\t\tA_SPECIES')
    suffix_out  = '/a_species.txt'
    filename_out= dataPath_out+continent+suffix_out
    # ISA-TAB COLUMNS: ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #                       |                  |                        |               |                        |           |                       |           |                                                          |                   |                                                                      |                |    
    headers     =           'Sample Name\t     Assay Name\t             Description\t   Protocol REF\t           Performer\t Protocol REF\t          Performer\t Characteristics [species assay result (VBcv:0000961)]\t    Term Source REF\t   Term Accession Number\t                                                Comment [ASSI]\n'
    values      = ''.join([ sampleName+'\t'+   sampleName+'.spid'+'\t'+ ''+'\t'+        speciesMethods1[i]+'\t'+ ''+'\t'+    speciesMethods2[i]+'\t'+''+'\t'+    species_raw_to_miro[speciesNames[i]]['sp_miro']+'\t'+      'MIRO'+'\t'+        species_raw_to_miro[speciesNames[i]]['id_miro'].split(':')[1]+'\t'+    assi[i]+'\n'     for i,sampleName in enumerate(sampleNames)]) # tab-delim text to parse
    #                       |                  |                        |               |                        |           |                       |           |                                                          |                   |                                                                      |                |
    #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    headers     = re.sub(r'\t ','\t',re.sub(r'( )\1{2,}', ' ', headers))
    parse       = headers+values
    file_out    = open(filename_out,'w')
    file_out.write(parse)
    file_out.close()

# DISPLAY UNIQUE METHODS        :       to help generate the i_investigations sheet
import pprint
pprint.pprint(('sample methods: ',[(i,j) for i,j in enumerate(unique_sample_methods)]))
pprint.pprint(('species methods: ',[(i,j) for i,j in enumerate(unique_species_methods)]))