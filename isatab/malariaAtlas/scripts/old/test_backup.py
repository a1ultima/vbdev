import numpy as np
import os 
import re 
import difflib

#===============#
# ONTOLOGY GEN  #
#===============#

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
        sp   = ' '.join(line_split[2].split())
        onto = ' '.join(line_split[1].split())
    except IndexError:
        continue
    species_to_ontoTerm[sp] = onto
file_in.close()
speciesNames_onto = species_to_ontoTerm.keys()

# COUNTRY_ID:ONTOTERMS         :       
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


#RENAME INCOMPATIBLE KEYS       :       e.g. spreadsheet has: "Anopheles (Cellia) flavirostris (Ludlow, 1914)" when the key is: "Anopheles flavirostris (Ludlow, 1914)"

# for species in species_to_ontoTerm.keys(): 
#     species_old = species
#     species_new = species.replace('Anopheles ','Anopheles (Cellia) ')
#     species_to_ontoTerm[species_new] = species_to_ontoTerm.pop(species_old)
# species_to_ontoTerm['Anopheles (Anopheles) barbirostris species complex']       = species_to_ontoTerm.pop('Anopheles (Cellia) barbirostris species complex')
# species_to_ontoTerm['Anopheles (Anopheles) sinensis species complex']           = species_to_ontoTerm.pop('Anopheles (Cellia) sinensis species complex')
# species_to_ontoTerm['Anopheles (Cellia) aconitus D\x9anitz, 1902']              = species_to_ontoTerm.pop('Anopheles (Cellia) aconitus Dönitz, 1902')
# species_to_ontoTerm['Anopheles (Anopheles) lesteri Baisas &amp; Hu, 1936 (formerly An. anthropophagus in China)'] = species_to_ontoTerm.pop('Anopheles (Cellia) lesteri Baisas &amp; Hu, 1936')
# species_to_ontoTerm['Anopheles (Cellia) merus D\xf6nitz, 1902']                 = species_to_ontoTerm.pop('Anopheles (Cellia) merus Dönitz, 1902')
# species_to_ontoTerm['Anopheles (Nyssorhynchus) darlingi Root, 1926']            = species_to_ontoTerm.pop('Anopheles (Cellia) darlingi Root, 1926')
# species_to_ontoTerm['Anopheles (Nyssorhynchus) nuneztovari species complex']    = species_to_ontoTerm.pop('Anopheles (Cellia) nuneztovari species complex')
# species_to_ontoTerm['Anopheles (Nyssorhynchus) albimanus Wiedemann, 1820']      = species_to_ontoTerm.pop('Anopheles (Cellia) albimanus Wiedemann, 1820')
# species_to_ontoTerm['Anopheles (Anopheles) pseudopunctipennis species complex'] = species_to_ontoTerm.pop('Anopheles (Cellia) pseudopunctipennis species complex')
# species_to_ontoTerm['Anopheles (Nyssorhynchus) albitarsis species complex']     = species_to_ontoTerm.pop('Anopheles (Cellia) albitarsis species complex')
# species_to_ontoTerm['Anopheles (Nyssorhynchus) marajoara Galv\xe3o &amp; Damasceno, 1942'] = species_to_ontoTerm.pop('Anopheles (Cellia) marajoara Galvão &amp; Damasceno, 1942')
# species_to_ontoTerm['Anopheles (Nyssorhynchus) aquasalis Curry, 1932']          = species_to_ontoTerm.pop('Anopheles (Cellia) aquasalis Curry, 1932')
# species_to_ontoTerm['Anopheles (Anopheles) quadrimaculatus Say, 1824']          = species_to_ontoTerm.pop('Anopheles (Cellia) quadrimaculatus Say, 1824')
# species_to_ontoTerm['Anopheles (Anopheles) freeborni Aitken, 1939']             = species_to_ontoTerm.pop('Anopheles (Cellia) freeborni Aitken, 1939')
# species_to_ontoTerm['Anopheles (Anopheles) sacharovi Favre, 1903']              = species_to_ontoTerm.pop('Anopheles (Cellia) sacharovi Favre, 1903')
# species_to_ontoTerm['Anopheles (Anopheles) atroparvus van Thiel, 1927']         = species_to_ontoTerm.pop('Anopheles (Cellia) atroparvus van Thiel, 1927')
# species_to_ontoTerm['Anopheles (Anopheles) messeae Falleroni, 1926']            = species_to_ontoTerm.pop('Anopheles (Cellia) messeae Falleroni, 1926')
# species_to_ontoTerm['Anopheles (Anopheles) labranchiae Falleroni, 1926']        = species_to_ontoTerm.pop('Anopheles (Cellia) labranchiae Falleroni, 1926')


#-------------------------------------------------------------------------------------------- 

#===============#
#  ISATAB GEN   #
#===============#

dataPath_out= '../data/isatab/'
dataPath_in = '../data/raw/'
suffix_in   = '_DVS.txt'
continents  = ['Asia','Africa','Americas','Europe']

print('generating isatabs...')
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
    data            = {}
    for col in lines_arr_T:
        data[col[0]]= col[1:]
    sampleNames     = list(data['id'])
    sampleMethods1  = list(data['sample method1'])
    sampleMethods2  = list(data['sample method2'])
    sampleMethods3  = list(data['sample method3'])
    sampleMethods4  = list(data['sample method4'])
    speciesMethods1 = list(data['id method1'])
    speciesMethods2 = list(data['id method2'])
    speciesNames    = [species.replace('"','') for species in list(data['species'])] #replace '"' with ''   :   some species are encapsulated in double-quotes
    countryIds      = list(data['country id'])
    latitudes       = list(data['latitude'])
    longitudes      = list(data['longitude'])
    year_starts     = list(data['year start'])
    year_ends       = list(data['year end'])
    month_starts    = list(data['month start'])
    month_ends      = list(data['month end'])
    dates           = [year+'-'+month_starts[i].zfill(2)+'-00'+'/'+year_ends[i]+'-'+month_ends[i].zfill(2)+'-00' for i,year in enumerate(year_starts)]


    citations       = list(data['citation'])
    assi            = list(data['ASSI'])

    # I_INVESTIGATIONS HELPER   :   finds all unique species ID & sampling methods, it accumulate across continents
    try:
        print('\t\tAccumulating methods: i_investigations...')
        unique_species_methods = list(set(list(data['id method1'])+list(data['id method1'])+unique_species_methods))
        unique_sample_methods  = list(set(list(data['sample method1'])+list(data['sample method2'])+list(data['sample method3'])+list(data['sample method4'])+unique_sample_methods))
    except NameError:
        print('\t\tAccumulating methods: i_investigations...')
        unique_species_methods = list(set(list(data['id method1'])+list(data['id method1'])))
        unique_sample_methods  = list(set(list(data['sample method1'])+list(data['sample method2'])+list(data['sample method3'])+list(data['sample method4'])))


# GENERATE ISATABS      :       s_samples, a_collections, a_species


#===============#
#   S_SAMPLES   #
#===============#
    print('\t\tGenerating s_samples sheet...')
    suffix_out  = '/s_samples.txt'
    filename_out= dataPath_out+continent+suffix_out

    # ISA-TAB COLUMNS: ----------------------------------------------------------------------------------------------------------------------------------
    #                           |                       |                   |               |               |                   |                       |
    headers     =               'Source Name\t          Sample Name\t       Description\t   Material Type\t Term Source REF\t   Term Accession Number\n'
    values      = ''.join([     'malariaAtlas'+'\t'+    sampleName+'\t'+    ''+'\t'+        ''+'\t'         ''+'\t'+            ''+'\n'                     for sampleName in sampleNames]) 
    #                           |                       |                    |              |               |                   |                       |
    #----------------------------------------------------------------------------------------------------------------------------------------------------

    headers     = re.sub(r'\t ','\t',re.sub(r'( )\1{2,}', ' ', headers))
    parse       = headers+values
    file_out = open(filename_out,'w')
    file_out.write(parse)
    file_out.close()

#===================#
#   A_COLLECTIONS   #
#===================#
    print('\t\tGenerating a_collections sheet...')
    suffix_out  = '/a_collections.txt'
    filename_out= dataPath_out+continent+suffix_out

    # ISA-TAB COLUMNS: --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #                           |                   |                               |               |                       |           |                       |           |                       |           |                       |           |               |                                                   |                   |                                                                   |                                                           |                                                   
    headers     =               'Sample Name\t      Assay Name\t                    Description\t   Protocol REF\t          Performer\t Protocol REF\t          Performer\t Protocol REF\t          Performer\t Protocol REF\t          Performer\t Date\t          Country\t                                           Term Source REF\t   Term Accession Number\t                                             Characteristics [Collection site latitude (VBcv:0000817)\t  Characteristics [Collection site longitude (VBcv:0000816)\n'
    values      = ''.join([     sampleName+'\t'+    sampleName+' collection'+'\t'+  ''+'\t'+        sampleMethods1[i]+'\t'+ ''+'\t'+    sampleMethods2[i]+'\t'+ ''+'\t'+    sampleMethods3[i]+'\t'+ ''+'\t'+    sampleMethods4[i]+'\t'+ ''+'\t'+    dates[i]+'\t'+  country_to_ontoTerm[countryIds[i]]['country']+'\t'+ 'GAZ'+'\t'+         country_to_ontoTerm[countryIds[i]]['gazTerm'].split(':')[1]+'\t'+   latitudes[i]+'\t'+                                          longitudes[i]+'\n'                                                  for i,sampleName in enumerate(sampleNames)]) # tab-delim text to parse
    #                           |                   |                               |               |                       |           |                       |           |                       |           |                       |           |               |                                                   |                   |                                                                   |                                                           |                                                   
    #--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #                                                                                                     ^--------------------^-----------------------^----------------------^ --> needs to be dictified 
    headers     = re.sub(r'\t ','\t',re.sub(r'( )\1{2,}', ' ', headers))
    parse       = headers+values
    file_out = open(filename_out,'w')
    file_out.write(parse)
    file_out.close()

#===============#
#   A_SPECIES   #
#===============#
    print('\t\tGenerating a_species sheet...')
    suffix_out  = '/a_species.txt'
    filename_out= dataPath_out+continent+suffix_out

    # ISA-TAB COLUMNS: --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #                        |                  |                           |               |                           |           |                       |           |                                                                       |                   |                                                                                                   |
    headers     =            'Sample Name\t     Assay Name\t                Description\t   Protocol REF\t              Performer\t Protocol REF\t          Performer\t Characteristics [species assay result (VBcv:0000961)]\t                 Term Source REF\t   Term Accession Number\n'
    values      = ''.join([  sampleName+'\t'+   sampleName+'.spid'+'\t'+    ''+'\t'+        speciesMethods1[i]+'\t'+    ''+'\t'+    speciesMethods2[i]+'\t'+''+'\t'+    difflib.get_close_matches(speciesNames[i],speciesNames_onto)[0]+'\t'+   'MIRO'+'\t'+        species_to_ontoTerm[difflib.get_close_matches(speciesNames[i],speciesNames_onto)[0]].split(':')[1]+'\n'       for i,sampleName in enumerate(sampleNames)]) # tab-delim text to parse
    #                        |                  |                           |               |                           |           |                       |           |                                                                       |                   |                                                                                                   |
    #--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    headers     = re.sub(r'\t ','\t',re.sub(r'( )\1{2,}', ' ', headers))
    parse       = headers+values
    file_out = open(filename_out,'w')
    file_out.write(parse)
    file_out.close()

# DISPLAY UNIQUE METHODS        :       to help generate the i_investigations sheet
import pprint
pprint.pprint(('sample methods: ',[(i,j) for i,j in enumerate(unique_sample_methods)]))
pprint.pprint(('species methods: ',[(i,j) for i,j in enumerate(unique_species_methods)]))