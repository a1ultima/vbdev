import numpy as np
import os 
import itertools
import difflib

print('Reading ontology dictionaries...')

dataPath_in = '../data/ontology/'
filename_in = dataPath_in+'MalariaAtlas_citations'
file_in     = open(filename_in,'r')

citation_to_pubmedId = {}
while True:
    line        = file_in.readline()
    line_split  = line.rstrip().split('\t') # Abdalla, H.M.A.B.   2008    18054056    Insecticide susceptibility and vector status of natural populations of Anopheles arabiensis from Sudan. Transactions of the Royal Society of Tropical Medicine and Hygiene, 102(3):263-71
    if line == "":                  # break when finished
        break        
    name    = line_split[0] # Abdalla, H.M.A.B.
    year    = line_split[1] # 2008
    pid     = line_split[2] # 18054056
    title   = line_split[3] # Insecticide susceptibility and vector status of natural populations of Anopheles arabiensis from Sudan.
    journal = line_split[4] # Transactions of the Royal Society of Tropical Medicine and Hygiene, 102(3):263-71
    citation_to_pubmedId[name+' '+year+' '+title+' '+journal] = pid #TODO: there are more unique citations in raw compared w/ Dan's ontologies
file_in.close()
citations_pubmed = citation_to_pubmedId.keys()

#--------------------------------------------------------

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
    while True: # WARNING: this does not skip @headers! 
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
    speciesNames = [species.replace('"','') for species in list(data[continent]['species'])] #replace '"' with ''   :   some species are encapsulated in double-quotes

#--------------------------------------------------------

# CITATIONS_RAW:(CITATIONS_PUBMED,ID_PUBMED)
print('Generating raw_citation:pubmed_citation dictionary...')
citations_raw_unique = list(set(list(itertools.chain(*[[citation.replace('"','') for citation in list(data[continent]['citation'])] for continent in continents]))))
citations_raw_unique = citations_raw_unique[1:] # this purges the header line away from raw (@headers!)
citation_raw_to_pubmed = {}
for count,citation_raw in enumerate(citations_raw_unique):
    print(count)
    try:
        citation_pubmed_match                  = difflib.get_close_matches(citation_raw, citations_pubmed, cutoff=0.4)[0]
        citation_raw_to_pubmed[citation_raw]   = {'sp_pubmed':citation_pubmed_match,'id_pubmed':citation_to_pubmedId[citation_pubmed_match]}
    except IndexError:
        print(citation_raw)

# # Takken, W., Geene, R., Adam, W., Jetten, T.H. and van der Velden, 
# # J.A. (2002).  <b>Distribution and dynamics of larval populations 
# # of <i>Anopheles messeae</i> and <i>A. atroparvus</i> in the delta
# # of the rivers Rhine and Meuse, The Netherlands.</b> <i>Ambio</i>, 
# # <b>31</b>(3):212-8
#Takken, W., Geene, R., Adam, W., Jetten, T.H. and van der Velden, J.A. (2002).  <b>Distribution and dynamics of larval populations of <i>Anopheles messeae</i> and <i>A. atroparvus</i> in the deltaof the rivers Rhine and Meuse, The Netherlands.</b> <i>Ambio</i>, <b>31</b>(3):212-8 
#--------------------------------------------------------






