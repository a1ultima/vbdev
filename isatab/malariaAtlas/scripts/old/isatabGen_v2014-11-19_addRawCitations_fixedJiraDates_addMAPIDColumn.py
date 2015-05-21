"""

Description: 

    Generate isatab sheets: s_samples.txt, a_collection.txt, a_species.txt per continent. Continents can be: Africa, Europe, Americas, Asia. 

    Essentially this involves two steps:
        (i) generating dictionaries that retrieve ontology terms given inputs from a single input spreadsheet.
        (ii) splitting the columns of a single input spreasheet (e.g. Europe_DVS.txt) into three separate sheets named above. 
        (iii) prints some data that may be useful for manually curating the i_investigations.txt sheet that completes the auto generated sheets above.

    data in: /home/ab108/0VB/isatab/malariaAtlas/data/raw/*_DVS.txt

    data out: /home/ab108/0VB/isatab/malariaAtlas/data/isatab/*

        where * is a continent e.g. Europe/**, where ** are three spreadsheets: a_collection.txt  a_species.txt  s_samples.txt that constitute the bulk of an isatab.

        Note: there is one more isatab that is required to complete the isatab sheets: i_investigation.txt which is manually curated. 

Pipeline Map:

    {{{ isatabGen*.py }}} --> isatabPrune*.py --> cp to /home/ab108/popbio/data_andy/isatabs/andy_0411_MAP/


    

"""

import pickle
import itertools
import numpy as np
import os 
import re 
import difflib
import re 

# Pipeline:
#    1. Generate dictionaries for: species, country, citations  --> species_to_miro, country_to_gaz, citations_to_pubmedId
#    2. Re-interpret the above dictionaries to allow conversion from raw_to_<aboveDictionary>  --> e.g. raw_to_species_to_miro
#    3. Read & Vectorize raw data columns
#    4. Parse vectorized columns partitioned into: s_samples, a_collection, a_species

# Functions:
def filterVector( vector ): 
    ''' Replaces elements: 'unknown' and 'not applicable' with ''
    '''
    out_vector = [x if ((not x=='unknown') and (not x=='not applicable')) else '' for x in vector]
    return out_vector

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
dataPath_onto   = '../data/ontology/'
suffixes_in     = ['country_id','species','MalariaAtlas_citations']

# SPECIES:ONTOTERMS         :  
print('\tSPECIESNAMES_to_MIROTERMS')     
filename_in = dataPath_onto+'species'   # /home/ab108/0VB/isatab/malariaAtlas/data/ontology/species
file_in     = open(filename_in)
species_to_ontoTerm = {}
while True:
    line = file_in.readline()
    if line == "":
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
filename_in = dataPath_onto+'country_id' 
file_in     = open(filename_in)  # /home/ab108/0VB/isatab/malariaAtlas/data/ontology/country_id
country_to_ontoTerm = {}
while True:
    line        = file_in.readline()
    line_split  = line.split('|')
    if line == "":
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

# CITATIONS:PUBMEDIDS
filename_in = dataPath_onto+'MalariaAtlas_citations' 
file_in     = open(filename_in,'r') # /home/ab108/0VB/isatab/malariaAtlas/data/ontology/MalariaAtlas_citations
citation_to_pubmedId = {}
while True:
    line        = file_in.readline()
    line_split  = line.rstrip().split('\t') # Abdalla, H.M.A.B.   2008    18054056    Insecticide susceptibility and vector status of natural populations of Anopheles arabiensis from Sudan. Transactions of the Royal Society of Tropical Medicine and Hygiene, 102(3):263-71
    if line == "":
        break    
    name    = line_split[0] # Abdalla, H.M.A.B.
    year1   = line_split[1] # 2008
    if '_' in line_split[2]:
        pid = ''
    else:
        pid = line_split[2] # 18054056
    title   = line_split[3] # Insecticide susceptibility and vector status of natural populations of Anopheles arabiensis from Sudan.
    journal = line_split[4] # Transactions of the Royal Society of Tropical Medicine and Hygiene, 102(3):263-71
    citation_to_pubmedId[name+' '+year1+' '+title+' '+journal] = pid #TODO: there are more unique citations in raw compared w/ Dan's ontologies
file_in.close()
citations_pubmed = citation_to_pubmedId.keys()

#-------------------------------------------------------------------------------------------- 

#===============#
#  ISATAB GEN   #
#===============#

dataPath_out= '../data/isatab/' # 
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
    #   READ DATA
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

    #   PAD MISSING DATA      :   Europe and Americas are missing columns in latter rows => uneven array => we pad rows to ensure even array by len(row) = 17
    
    print('\t\tPadding missing data...')

    lines_even = []
    for row in lines:
        if len(row) < 17:
            row = row+['']*(17-len(row))
        lines_even.append(row)

    #   Convert lines of data into vectors of n elements
    print('\t\tVectorizing data...')

    lines_arr       = np.array(lines_even)
    lines_arr_T     = lines_arr.T           # rows of data that are read in are transposed because we want to deal with the data in terms of columns for later

    data[continent] = {}

    for col in lines_arr_T:                 # col[0]  = headers // e.g. 'latitude' // see: /home/ab108/0VB/isatab/malariaAtlas/data/raw/Europe_DVS.txt
        data[continent][col[0]]= col[1:]    # col[1:] = values  // e.g. 37.75

# RE-INTERPRET ONTOLOGIES:

#       SPECIES_RAW:(SPECIES_MIRO,ID_MIRO)

try: 
    species_raw_to_miro = pickle.load(open(dataPath_onto+'species_raw_to_miro.p','r'))
    print('Importing raw_species:miro_species dictionary...')
except IOError:
    print('Generating raw_species:miro_species dictionary...')
    species_raw_unique  = list(set(list(itertools.chain(*[[species.replace('"','') for species in list(data[continent]['species'])] for continent in continents]))))
    species_raw_to_miro = {}
    for species_raw in species_raw_unique:
        species_miro_match                  = difflib.get_close_matches(species_raw,sp_miro)[0]
        species_raw_to_miro[species_raw]    = {'sp_miro':species_miro_match,'id_miro':species_to_ontoTerm[species_miro_match]}
    pickle.dump(species_raw_to_miro,open(dataPath_onto+'species_raw_to_miro.p','w'))

#       CITATIONS_RAW:(CITATIONS_PUBMED,ID_PUBMED)

try:
    citation_raw_to_pubmed = pickle.load(open(dataPath_onto+'citation_raw_to_pubmed.p','r'))
    print('Importing raw_citation:pubmed_citation dictionary...')
except IOError:
    print('Generating raw_citation:pubmed_citation dictionary...')
    citations_raw_unique = list(set(list(itertools.chain(*[[citation.replace('"','') for citation in list(data[continent]['citation'])] for continent in continents]))))
    citations_raw_unique = citations_raw_unique[1:] # this purges the header line away from raw (@headers!)
    citation_raw_to_pubmed = {}
    for count,citation_raw in enumerate(citations_raw_unique):
        print count
        citation_pubmed_match                  = difflib.get_close_matches(citation_raw, citations_pubmed, cutoff=0.4)[0]
        citation_raw_to_pubmed[citation_raw]   = {'citation_pubmed':citation_pubmed_match,'id_pubmed':citation_to_pubmedId[citation_pubmed_match]}
    citation_raw_to_pubmed[''] = {'id_pubmed':'Unknown','citation_pubmed':'Unknown'} # WARNING: append a key for "unknown" where citations are not available
    pickle.dump(citation_raw_to_pubmed,open(dataPath_onto+'citation_raw_to_pubmed.p','w'))


# WRITE ISATABS      :       s_samples, a_collection, a_species

unique_species_methods  = [] # list for gathering names of unique species/sampling methods, these become useful for generating the i_investigations.txt sheet for each continent: e.g. /home/ab108/popbio/data_andy/isatabs/andy_0411_MAP/old/Europe/i_investigation.txt
unique_sample_methods   = []


print('writing ISA-TABs to file...')
for continent in continents:
    print('\t'+continent)

    ######################
    # FORMAT DATA COLUMNS:       each line here is a columns of data that will later be split between: s_samples, a_collection, a_species
    ######################

    sampleNames     = list(data[continent]['id'])
    sampleMethods1  = filterVector(list(data[continent]['sample method1']))
    sampleMethods2  = filterVector(list(data[continent]['sample method2']))
    sampleMethods3  = filterVector(list(data[continent]['sample method3']))
    sampleMethods4  = filterVector(list(data[continent]['sample method4']))
    speciesMethods1 = filterVector(list(data[continent]['id method1']))
    speciesMethods2 = filterVector(list(data[continent]['id method2']))

    speciesNames    = [species.replace('"','') for species in list(data[continent]['species'])]     #replace '"' with ''   :   some species are encapsulated in double-quotes
    countryIds      = list(data[continent]['country id'])
    latitudes       = list(data[continent]['latitude'])
    longitudes      = list(data[continent]['longitude'])

    # Dates:
    year_starts     = list(data[continent]['year start'])
    year_ends       = list(data[continent]['year end'])
    month_starts    = list(data[continent]['month start'])
    month_ends      = list(data[continent]['month end'])
    dates           = [year+'-'+month_starts[i].zfill(2)+'/'+year_ends[i]+'-'+month_ends[i].zfill(2) for i,year in enumerate(year_starts)]
    
    #dates           = [i.replace('-00/-00','') for i in dates] # @DATEBUG: there are some date that have missing years but still have their months, e.g. -06/-09. We need to get rid of them just as we did with empty -00/-00 dates 
    dates           = [re.sub('-[0-9][0-9]/-[0-9][0-9]','',i) for i in dates ] # make a regex replacement, e.g. -[0-9][0-9]/-[0-9][0-9]

    #citations       = [citation.replace('"','') for citation in list(data[continent]['citation'])]  #replace '"' with ''
    citations       = list(data[continent]['citation'])
    assi            = list(data[continent]['ASSI'])

    ##############################
    # I_INVESTIGATIONS HELPER   :  I will use this to help me annotate the i_investigations.txt sheet: the following code finds all unique species ID & sampling methods, it accumulate across continents
    ##############################
    

    try:
        print('\t\tAccumulating methods: i_investigations...OK')
        print(data[continent].keys())

        unique_species_methods = list(set(list(data[continent]['id method1'])+list(data[continent]['id method2'])+unique_species_methods))
        unique_sample_methods  = list(set(list(data[continent]['sample method1'])+list(data[continent]['sample method2'])+list(data[continent]['sample method3'])+list(data[continent]['sample method4'])+unique_sample_methods))
    
    except NameError as detail:
        print('\t\tSOMETHING IS WRONG!!!')

        print detail

        unique_species_methods = list(set(list(data[continent]['id method1'])+list(data[continent]['id method2'])))
        unique_sample_methods  = list(set(list(data[continent]['sample method1'])+list(data[continent]['sample method2'])+list(data[continent]['sample method3'])+list(data[continent]['sample method4'])))




#############################################
# WRITE OUT DATA TO 3 ISATAB.txt SHEETS: (i) s_samples, (ii) 
#############################################


    #=======================#
    #  s_samples.txt        #
    #=======================#
    print('\t\tGenerating Sheets...')
    print('\t\t\tS_SAMPLES')
    suffix_out  = '/s_samples.txt'
    filename_out= dataPath_out+continent+suffix_out
    # ISA-TAB COLUMNS: ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #                       |                       |                       |               |                   |                   |                       |                                   |                                                                        |                                    
    headers     =           'Source Name\t          Sample Name\t           Description\t   Material Type\t     Term Source REF\t   Term Accession Number\t Comment [citation]\t                 Comment [PubMed ID]\t                                                   Comment [MAP ID]\n'                  
    values      = ''.join([ 'malariaAtlas'+'\t'+    'MAP.'+sampleName+'\t'+ ''+'\t'+        'population'+'\t'   'OBI'+'\t'+         '0000181'+'\t'+         citations[i].replace('"','')+'\t'+   citation_raw_to_pubmed[citations[i].replace('"','')]['id_pubmed']+'\t'+ 'MAP.'+sampleName+'\n'                  for i,sampleName in enumerate(sampleNames)]) 
    #                       |                       |                       |               |                   |                   |                       |                                   |                                                                        |                                    
    #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    headers     = re.sub(r'\t ','\t',re.sub(r'( )\1{2,}', ' ', headers))

    assert len(headers)==2, "there are multiple newline characters in headers, CTRL+F: 'headers     = '"

    parse       = headers+values
    file_out    = open(filename_out,'w')
    file_out.write(parse)
    file_out.close()

    #=======================#
    #   a_collection.txt    #
    #=======================#
    print('\t\t\tA_collection')
    suffix_out  = '/a_collection.txt'
    filename_out= dataPath_out+continent+suffix_out
    # ISA-TAB COLUMNS: ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #                                   |                       |                                       |               |                       |                       |                       |                       |               |                                                       |                   |                                                                   |                                                               |                                                           |
    headers     = re.sub(r' {2,}','',   'Sample Name\t          Assay Name\t                            Description\t   Protocol REF\t          Protocol REF\t          Protocol REF\t          Protocol REF\t          Date\t          Characteristics [Collection site (VBcv:0000831)]\t      Term Source REF\t   Term Accession Number\t                                             Characteristics [Collection site latitude (VBcv:0000817)]\t     Characteristics [Collection site longitude (VBcv:0000816)]\n')
    values      = ''.join([             'MAP.'+sampleName+'\t'+ 'MAP.'+sampleName+'.collection'+'\t'+   ''+'\t'+        sampleMethods1[i]+'\t'+ sampleMethods2[i]+'\t'+ sampleMethods3[i]+'\t'+ sampleMethods4[i]+'\t'+ dates[i]+'\t'+  country_to_ontoTerm[countryIds[i]]['country']+'\t'+     'GAZ'+'\t'+         country_to_ontoTerm[countryIds[i]]['gazTerm'].split(':')[1]+'\t'+   latitudes[i]+'\t'+                                              longitudes[i]+'\n'                                           for i,sampleName in enumerate(sampleNames)]) # tab-delim text to parse
    #                                   |                       |                                       |               |                       |                       |                       |                       |               |                                                       |                   |                                                                   |                                                               |                                                           |
    #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #                                                                                                     ^--------------------^-----------------------^----------------------^ --> needs to be dictified 
    headers     = re.sub(r'\t ','\t',re.sub(r'( )\1{2,}', ' ', headers))
    parse       = headers+values
    file_out    = open(filename_out,'w')
    file_out.write(parse)
    file_out.close()

    #======================#
    #   a_species.txt      #
    #======================#
    print('\t\t\tA_SPECIES')
    suffix_out  = '/a_species.txt'
    filename_out= dataPath_out+continent+suffix_out
    # ISA-TAB COLUMNS: ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #                       |                       |                               |               |                        |                       |                                                          |                   |                                                                      |                |    
    headers     =           'Sample Name\t          Assay Name\t                    Description\t   Protocol REF\t           Protocol REF\t          Characteristics [species assay result (VBcv:0000961)]\t    Term Source REF\t   Term Accession Number\t                                                Comment [ASSI]\n'
    values      = ''.join([ 'MAP.'+sampleName+'\t'+ 'MAP.'+sampleName+'.spid'+'\t'+ ''+'\t'+        speciesMethods1[i]+'\t'+ speciesMethods2[i]+'\t'+species_raw_to_miro[speciesNames[i]]['sp_miro']+'\t'+      'MIRO'+'\t'+        species_raw_to_miro[speciesNames[i]]['id_miro'].split(':')[1]+'\t'+    assi[i]+'\n'     for i,sampleName in enumerate(sampleNames)]) # tab-delim text to parse
    #                       |                       |                               |               |                        |                       |                                                          |                   |                                                                      |                |
    #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    headers     = re.sub(r'\t ','\t',re.sub(r'( )\1{2,}', ' ', headers))
    parse       = headers+values
    file_out    = open(filename_out,'w')
    file_out.write(parse)
    file_out.close()



# DISPLAY UNIQUE METHODS        :       to help generate the i_investigations sheet
import pprint
pprint.pprint(('sample methods: ',[(i,j) for i,j in enumerate(unique_sample_methods)]))
pprint.pprint(('species methods: ',[(i,j) for i,j in enumerate(unique_species_methods)]))