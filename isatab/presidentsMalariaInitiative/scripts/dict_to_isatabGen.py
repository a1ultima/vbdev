# coding=utf-8


#pdb.set_trace() #@1212
###########
# Imports #
###########

from os import listdir
import pickle
import numpy as np
import os 
import difflib

import csv 
import numpy as np
import pdb
from copy import copy

#############
# Functions #
#############

def check_path_exists_else_make_it( path ):
    """ Check a path exists, if not make it. 

    Returns:
        None
    Args:
        path:   relative path. e.g. dir = '..\Documents\Other'
    Exceptions:
           None

    """
    import shutil

    if os.path.exists(path):
        pass
    else:
        os.makedirs(path)

def filter_assay_rows( assay_data_rowise, assay_string_match ):

    """ remove all rows of data not pertaining to the desired assay 

    Args:
        assay_string_match:     e.g. "WHO test kit_adults"  (from: header_to_datacolumn["Test                                           type"]["raw_dataset_column"])
        assay_data_rowise:      e.g. a_IR_WHO ( = a_IR_WHO.T ): e.g. ['PMI.96', 'PMI.IR_WHO.96insecticide', 'IR_WHO', '', '', '', 'lambda-cyhalothrin', 'MIRO', '10000125', '0.05', 'percent', 'UO', '0000187', '0000032', '0000032', 'UO', '0000032', 'Raw Data File']
    Returns:
        similar to assay_data_rowise, but with rows only matching to those relevant to the assay specified in assay_string_match

    """
    # get all row pertaining to IR_WHO assays
    assay_specification_column = header_to_datacolumn["Test                                           type"]["raw_dataset_column"]

    # get indices of rows of data not relevant to IR_WHO assays
    row_indices_matching_specified_assay    = [i for i,row in enumerate(assay_specification_column) if not row==assay_string_match]

    filtered_assay_data                     = np.delete(assay_data_rowise,row_indices_matching_specified_assay,axis=0)

    return filtered_assay_data

########
# Main #
########

#
# Check data paths exist
#

# check if necessary paths exist, else make them
#
# @DONE: check if file exists else create it
check_path_exists_else_make_it("../data/")
check_path_exists_else_make_it("../data/ontologies")
check_path_exists_else_make_it("../data/raw")
check_path_exists_else_make_it("../data/isatab")

#
# ...
#

import time
start_time = time.time()

#
# Read .csv file representing the punett
#
reader      = csv.reader(open("../data/raw/IR Database_PMI Dataset 22122014.csv","rb"),delimiter=",")
x           = list(reader)
dataset     = np.array(x)
superheaders= dataset[0]
headers     = dataset[1] 

dataset_transposed = list(dataset[2:].T)  # so that columns become stacked as rows

#
# 1. Get some ontologies 
#
#       search for all the pickle files (which should be ontology keys) 
#
ontologyKey_to_variableName = {}

pickle_names = [fn for fn in listdir("../data/ontologies/") if ".p" in fn]

for pickle_name in pickle_names:
    with open('../data/ontologies/'+pickle_name, 'rb') as handle:
        pickled_ontology = pickle.load(handle)
        pickle_name = pickle_name.replace(".p","")
        ontologyKey_to_variableName[pickle_name] = pickled_ontology  # e.g. ../data/ontologies/country_to_gaz.txt

#
# 2. Split the columns into sheets they belong to 
#
#       s_samples
#
# put the transposed data into a dict

header_to_datacolumn = {}

for i,header in enumerate(headers):
    header_to_datacolumn[header] = { 'raw_dataset_column':dataset_transposed[i],\
                                     'mapped_ontology_columns':{} }

# test it DONE

# for a given data_transposed column, get all the mapped onto terms in the ontologies dict, as their own columns

# get the ontologies out one by one, then get their corresponding columns of data out

# create a dictionary that translates between keys in header_to_datacolumn vs. keys in ontologyKey_to_variableName
headerToDatacolumnKEY_to_ontologyKeyToVariableNameKEY = {
    'Species tested':       'species_to_miroTerm_miroId',
    'Insecticide tested':   'insecticide_to_ontoTerms',
    'Country':              'country_to_gazTerm_gazId',
    'Start month':          'startMonth_to_numMonth',
    'End month':            'endMonth_to_numMonth',
    'Year':                 'yearStartToEnd_yearStart_yearEnd',
    'Stage tested (and origin)': 'stageTested_ontoTermAdult_termSourceRef_accnNum'
}


#
# Get onto term columns 
#   For every column with a matching ontology, parse the raw values through dictionaries to obtain corresponding ontology values/terms
#
for key1 in headerToDatacolumnKEY_to_ontologyKeyToVariableNameKEY.keys():

    # e.g. key1: 'Species tested', e.g. key2: 'species_to_miroTerm_miroId' 
    #   1. species 
    key1 = key1 # e.g. 'Species tested', for more see above: CTRL+F: "e.g. key1:"
    key2 = headerToDatacolumnKEY_to_ontologyKeyToVariableNameKEY[key1]  # e.g. key2: 'species_to_miroTerm_miroId' 

    for row in header_to_datacolumn[key1]['raw_dataset_column']:

        if "\xe6g" in row:
            row = row.replace("\xe6g","\xc3\xa6g")

        try: 
            # IF there is a second layer of keys for the dict, then iterate over those as subkeys
            subkeys = ontologyKey_to_variableName[key2][row].keys()
        except AttributeError:
            # ELSE there is only one layer of keys for the dict, then just treat that as a pointless layer of keys to make the dict symmetrical
            threeKeyLevelClone = copy(ontologyKey_to_variableName[key2])

            for toBeLowerKey in threeKeyLevelClone:
                keyClone = copy(threeKeyLevelClone[toBeLowerKey])
                ontologyKey_to_variableName[key2][toBeLowerKey] = {'justForSymmetry':keyClone}

        # @SUBKEY:  ...@todo...
        for ontology_term in ontologyKey_to_variableName[key2][row].keys():

            # Grab the ontology term's value that is mapped to "<row>" ...
            ontology_value = ontologyKey_to_variableName[key2][row][ontology_term] # e.g. Ghana, GAZ, 1000010

            try:
                header_to_datacolumn[key1]['mapped_ontology_columns'][ontology_term].append(ontology_value)
            except KeyError:
                header_to_datacolumn[key1]['mapped_ontology_columns'][ontology_term] = [ontology_value]                


############################################################################
############################################################################
# Generate ISA-Tab's tabs: s_samples.txt, a_collections.txt, g_assay...    #
############################################################################
############################################################################

###################
## s_samples.txt ## @s_samples
###################

# @todo: test s_samples for errors then show Bob
sample_names                = header_to_datacolumn["ORIG SOURCE REF"]["raw_dataset_column"]

# 2. Create a source name column    header: "Source Name"
col_sample_names            = sample_names

nrows = len(col_sample_names)
col_source_names            = np.array(["PresidentsMalariaInitiative"]*nrows)

# 1. Create sample name column:     header: "Sample Name"  foreign-key:one-to-one:@s_samples:a_collections, see: "@a_collections:s_samples"

# 3. header: "Description"
col_description             = np.array([""]*nrows)

# X. header: "Comment [ORIG REF]"
col_comment_origRef         = np.array(header_to_datacolumn["ORIG REF#"]["raw_dataset_column"])

# ...
col_comment_ref             = np.array(header_to_datacolumn["REF#"]["raw_dataset_column"])

# ...
col_comment_origSourceRef   = np.array(header_to_datacolumn["ORIG SOURCE REF"]["raw_dataset_column"])

# ...
col_comment_ref             = np.array(["PMI SOURCE REF#"]*nrows)

col_comment                 = np.array([""]*nrows)

# 4. header: "Material Type"     
col_material_type           = np.array(["pool"]*nrows)  # according to ISA-Tab:fonseca:kathrin: ISA-Tab: https://docs.google.com/spreadsheets/d/16_UQVMct4b6p3KmQe2-jftBnCgrvLTS25XwxaVvn5Qw/edit#gid=1341677884

# 5. header: "Term Source Ref"
col_term_source_ref         = np.array(["OBI"]*nrows)

# 6. header: "Term Accession Number"
col_term_accession_number   = np.array(["0000181"]*nrows)

# 10. header: "Comment [age]"
#col_age                     = np.array([""]*nrows)

# 12. header: "Characteristics [sex (EFO:0000695)]"
col_sex                     = np.array(["female"]*nrows)

# 13. header: "Term Source Ref"
col_sex_termSourceRef       = np.array(["PATO"]*nrows)

# 14. header: "Term Accession Number"
col_sex_accn                = np.array(["0000383"]*nrows)

# 7. header: 'Comment [Stage tested (and origin)]'
col_stage_tested            = np.array(header_to_datacolumn['Stage tested (and origin)']['raw_dataset_column'])

# 8. header: 'Characteristics [developmental stage (EFO:0000399)]'
col_stage_ontoTerm          = np.array(header_to_datacolumn['Stage tested (and origin)']['mapped_ontology_columns']['ontoTermAdult'])

# 9. header: 'Term Source Ref' 
col_stage_termSourceRef     = np.array(header_to_datacolumn['Stage tested (and origin)']['mapped_ontology_columns']['termSourceRef'])

# 10. header: 'Term Accession Number'
col_stage_accn              = np.array(["0000655"]*nrows)

#pdb.set_trace() # @1552

# 10. header: 'Characteristics [age (EFO:0000246]'
col_ageEFO                  = np.array(["2-5"]*nrows)

# 10. header: 'Unit'
col_ageEFO_unit             = np.array(["day"]*nrows)

# 10. header: 'Term Source Ref'
col_ageEFO_termSourceRef    = np.array(["UO"]*nrows)

# 10. header: 'Term Accession Number'
col_ageEFO_termAccnNo       = np.array(["0000033"]*nrows)

# @unite-columns
s_samples_headers           = np.array([    'Source Name',\
                                            'Sample Name',\
                                            'Description',\
                                            "Comment [PMI ORIG REF#]",\
                                            "Comment [PMI REF#]",\
                                            "Comment [PMI SOURCE REF#]",\
                                            'Comment [comment]',\
                                            'Material Type',\
                                            'Term Source Ref',\
                                            'Term Accession Number',\
                                            'Characteristics [sex (EFO:0000695)]',\
                                            'Term Source Ref',\
                                            'Term Accession Number',\
                                            'Comment [Stage tested (and origin)]',\
                                            'Characteristics [developmental stage (EFO:0000399)]',\
                                            'Term Source Ref',\
                                            'Term Accession Number',\
                                            'Characteristics [age (EFO:0000246)]',\
                                            'Unit',\
                                            'Term Source Ref',\
                                            'Term Accession Number' ])

### stack the s_samples columns vertically as rows 
#   [col 1: ..... ]
#   [col 2: ..... ]
#   [ ............]

s_samples = np.vstack(( col_source_names,\
                        col_sample_names,\
                        col_description,\
                        col_comment_origRef,\
                        col_comment_ref,\
                        col_comment_origSourceRef,\
                        col_comment,\
                        col_material_type,\
                        col_term_source_ref,\
                        col_term_accession_number,\
                        col_sex,\
                        col_sex_termSourceRef,\
                        col_sex_accn,\
                        col_stage_tested,\
                        col_stage_ontoTerm,\
                        col_stage_termSourceRef,\
                        col_stage_accn,\
                        col_ageEFO,\
                        col_ageEFO_unit,\
                        col_ageEFO_termSourceRef,\
                        col_ageEFO_termAccnNo ))

### transpose such that columns fall into columns
# [ col1 col2 ... ]
# [ ..............]
# [ ..............]
s_samples = s_samples.T

### stack the column headers on top
s_samples = np.vstack((s_samples_headers,s_samples))


# @save:s_samples
np.savetxt("../data/isatab/s_samples.tsv",      s_samples,      delimiter="\t", fmt="%s")


#######################
## a_collections.txt ## @a_collections
#######################

# 1. look at kathrin / fonseca IR study ISA-Tab and get all column headers below \/
## the following are column headers of @Kathrin's ISA-Tab (@fonseca): a_collections.txt

# Sample Name
col_collection_sample_names           = np.array([i for i in sample_names])  # foreign-key:one-to-one:see:"s_samples:a_collections"

# Assay Name
col_collection_assay_names            = np.array(["PMI.collection."+i for i in sample_names])  # we add a new "PMI.collection.i" suffix to map, assay names, one-to-one, with PMI.i sample names @foreign-key:one-to-one:@a_collections:s_samples see: "@s_samples:a_collections"

# Description
### e.g.: "human landing patch"
col_collection_description            = np.array(["catch of live specimens"]*nrows)

# Protocol REF
### e.g.: "CATCH"
col_collection_protocolRef            = np.array(["CATCH"]*nrows)

# Performer
### e.g.: " "
col_collection_performer              = np.array([" "]*nrows)

# Date
### e.g.: "2006-05/2006-11"
col_collection_date                   = np.array([" "]*nrows)
months_start                          = header_to_datacolumn['Start month']['raw_dataset_column']  # @raw-dataset @main @dataset
months_end                            = header_to_datacolumn['End month']['raw_dataset_column']

# obtain two columns: year_start, year_end from one column "year": [ (x), (y,y), ... ] --> [ (x,x), (y,y), ... ]
years_start = []  # e.g. ...
years_end   = []  # e.g. [ (x), (y,y), ... ] --> [ (x,x), (y,y), ... ]
year_rows   = header_to_datacolumn['Year']['raw_dataset_column']

for year_row in year_rows:
    if "-" in year_row:
        # sometimes year_row looks like this: 2001-2002, so we split them:
        year_start_end = year_row.split("-")
        years_start.append(year_start_end[0])
        years_end.append(year_start_end[1])

    elif "NR" in year_row:
        # but sometimes there is no available date, so we put an empty string there
        years_start = ""
        years_end = ""

    else: 
        # other times they look like this: 2001, so we copy them:
        years_start.append(year_row)
        years_end.append(year_row)

# dict mapping string month to numerical month
strMonth_to_numMonth = {
    "Jan":"01",
    "Feb":"02",
    "Mar":"03",
    "Apr":"04",
    "May":"05",
    "Jun":"06",
    "Jul":"07",
    "Aug":"08",
    "Sep":"09",
    "Oct":"10",
    "Nov":"11",
    "Dec":"12",
    "NR":"NR",
    "Jul, Aug":"07"
}

# generate some dates in the format required by ISA-Tab 
col_collection_date = []

for i,month_start in enumerate(months_start):
    month_start = strMonth_to_numMonth[months_start[i]]
    month_end   = strMonth_to_numMonth[months_end[i]]
    year_start  = years_start[i]
    year_end    = years_end[i]
    date        = year_start+"-"+month_start+"/"+year_end+"-"+month_end
    date        = date.replace("-NR","")  # sometimes the date is not included, represented as "-NR", these need to be trimmed off
    col_collection_date.append(date)

# save it as a column 
col_collection_date             = col_collection_date

# Characteristics [sampling time (EFO:0000689)]
### e.g.: "17:30/20:00"
# @DONE: check sampling time from the report in the email > find email from WL
col_collection_sampling_time    = [""]*nrows  # left empty since we do not know sampling time

# Characteristics [Temperature at time of collection (EFO:0001702)]
### e.g.: "25/28"
col_collection_temperature      = [""]*nrows  # left empty since temperature unavailable

# Unit
## e.g.: "degree Celsius"
col_collection_temperature_unit = [""]*nrows  # left empty since temperature unavailable

# Term Source Ref
## e.g.: "UO"
col_collection_temperature_source_ref = [""]*nrows  # left empty since temperature unavailable

# Term Accession Number
## e.g.: "0000027"
col_collection_temperature_accn = [""]*nrows  # left empty since temperature unavailable

# Comment [household ID]
## e.g.: " "
col_collection_householdId      = [""]*nrows  # left empty since household ID unavailable

# Characteristics [Collection site (VBcv:0000831)]
## e.g.: "Municipality of Antioquia"
col_collection_site             = np.array(header_to_datacolumn['Country']['mapped_ontology_columns']['gazTerm'])

# Term Source Ref
## e.g.: "GAZ"
col_collection_site_termSourceRef = np.array(header_to_datacolumn['Country']['mapped_ontology_columns']['source'])

# Term Accession Number
## e.g.: "00073117"
col_collection_site_accn          = np.array(header_to_datacolumn['Country']['mapped_ontology_columns']['gazId'])

# Characteristics [Collection site latitude (VBcv:0000817)]
## e.g.: "7.582777778"
col_collection_latitude         = header_to_datacolumn['Latitude UTM_X UPDATE2']['raw_dataset_column']

# Characteristics [Collection site longitude (VBcv:0000816)]
## e.g.: "-75.35194444"
col_collection_longitude        = header_to_datacolumn['Longitude UTM_Y UPDATE2']['raw_dataset_column']

# Characteristics [Collection site altitude (VBcv:0000832)]
## e.g.: ""
col_collection_altitude         = [""]*nrows  # not available

# Characteristics [Collection site location (VBcv:0000698)]
## e.g.: ""
col_collection_location         = [""]*nrows  # not available

# Characteristics [Collection site village (VBcv:0000829)]
## e.g.: ""
col_collection_village          = [""]*nrows  # not available

# Characteristics [Collection site locality (VBcv:0000697)]
## e.g.: ""
col_collection_locality         = np.array(header_to_datacolumn["Commune 3rd admin level"]["raw_dataset_column"])

# Characteristics [Collection site suburb (VBcv:0000845)]
## e.g.: ""
col_collection_suburb           = np.array([""]*nrows)

# Characteristics [Collection site city (VBcv:0000844)]
## e.g.: ""
col_collection_city             = np.array([""]*nrows)

# Characteristics [Collection site county (VBcv:0000828)]
## e.g.: ""
col_collection_county           = np.array([""]*nrows)

# Characteristics [Collection site district (VBcv:0000699)]
# @DONE: consolidate 'District Update' w/ 'District 2nd admin level' by prioritising the latter if available
districts              = header_to_datacolumn['District Update']['raw_dataset_column']
districts_updated      = header_to_datacolumn['District 2nd admin level']['raw_dataset_column']
districts_consolidated = []
for i,district in enumerate(districts):
    district_updated = districts_updated[i]
    if districts_updated[i]=="":
        districts_consolidated.append(district)
    else:
        districts_consolidated.append(district_updated)
col_collection_district = districts_consolidated

# Characteristics [Collection site province (VBcv:0000700)]
countries = header_to_datacolumn['Country']['mapped_ontology_columns']['gazTerm']
provinces = header_to_datacolumn['Province             1st admin level']['raw_dataset_column']
non_provinces = ["Southern, Eastern, Northern, Western"]

# since some provinces are stupidly named "Southern", the following consolidates the Country and Province columns into "Southern Angola". 
col_collection_provinces = []
for i,province in enumerate(provinces):
    if province in non_provinces:
        province = province+" "+countries[i]
    else:
        col_collection_provinces.append(province)
# test if something went horribly wrong above
assert len(col_collection_provinces)==len(provinces)==nrows 

# Characteristics [Collection site country (VBcv:0000701)]
col_collection_country = countries

# a_collection headers (from: Fonseca ISA-Tab: https://goo.gl/607nLN)
a_collection_headers = [    "Sample Name",\
                            "Assay Name",\
                            "Description",\
                            "Protocol REF",\
                            "Performer",\
                            "Date",\
                            "Characteristics [sampling time (EFO:0000689)]",\
                            "Characteristics [Temperature at time of collection (EFO:0001702)]",\
                            "Unit",\
                            "Term Source Ref",\
                            "Term Accession Number",\
                            "Comment [household ID]",\
                            "Characteristics [Collection site (VBcv:0000831)]",\
                            "Term Source Ref",\
                            "Term Accession Number",\
                            "Characteristics [Collection site latitude (VBcv:0000817)]",\
                            "Characteristics [Collection site longitude (VBcv:0000816)]",\
                            "Characteristics [Collection site altitude (VBcv:0000832)]",\
                            "Characteristics [Collection site location (VBcv:0000698)]",\
                            "Characteristics [Collection site village (VBcv:0000829)]",\
                            "Characteristics [Collection site locality (VBcv:0000697)]",\
                            "Characteristics [Collection site suburb (VBcv:0000845)]",\
                            "Characteristics [Collection site city (VBcv:0000844)]",\
                            "Characteristics [Collection site county (VBcv:0000828)]",\
                            "Characteristics [Collection site district (VBcv:0000699)]",\
                            "Characteristics [Collection site province (VBcv:0000700)]",\
                            "Characteristics [Collection site country (VBcv:0000701)]"]
#
# @unite-columns:a_species
#
a_collection = np.vstack((  col_collection_sample_names,\
                            col_collection_assay_names,\
                            col_collection_description,\
                            col_collection_protocolRef,\
                            col_collection_performer,\
                            col_collection_date,\
                            col_collection_sampling_time,\
                            col_collection_temperature,\
                            col_collection_temperature_unit,\
                            col_collection_temperature_source_ref,\
                            col_collection_temperature_accn,\
                            col_collection_householdId,\
                            col_collection_site,\
                            col_collection_site_termSourceRef,\
                            col_collection_site_accn,\
                            col_collection_latitude,\
                            col_collection_longitude,\
                            col_collection_altitude,\
                            col_collection_location,\
                            col_collection_village,\
                            col_collection_locality,\
                            col_collection_suburb,\
                            col_collection_city,\
                            col_collection_county,\
                            col_collection_district,\
                            col_collection_provinces,\
                            col_collection_country  ))
#
# rows are flipped to columns
#
a_collection = a_collection.T

### stack the column headers on top
a_collection = np.vstack((a_collection_headers,a_collection))

# @save:a_collections
np.savetxt("../data/isatab/a_collection.tsv",   a_collection,   delimiter="\t", fmt="%s")


#############
# a_species #  @@A_sPECIES
#############

#Sample Name
# "Corrales-El Playón"
col_species_sample_names        = np.array([i for i in sample_names])

#Assay Name
# "Corrales-El"
col_species_assay_names         = np.array(["PMI.species."+i for i in sample_names])

#Description
# "taxonomical identification"
col_species_description         = np.array([""]*nrows)

#Protocol REF
# "SPECIES"
col_species_protocolRef         = np.array(["SPECIES"]*nrows)

#Protocol REF
# "SPECIES"
col_species_protocolRef         = np.array(["SPECIES_PCR"]*nrows)

#Performer
# ""
col_species_performer           = np.array([""]*nrows)

#Date
# ""
col_species_date                = np.array([""]*nrows)  # @check:Fonseca vs. Uniqeu-PMI

#Characteristics [species assay result (VBcv:0000961)]
# "Anopheles nuneztovari"
col_species_result_value        = np.array(header_to_datacolumn["Species tested"]["mapped_ontology_columns"]["miroTerm"])

#Term Source Ref
# "MIRO"
col_species_result_termSourceRef= np.array(["MIRO"]*nrows)

#Term Accession Number
# "40000165"
col_species_result_accn         = np.array(header_to_datacolumn["Species tested"]["mapped_ontology_columns"]["miroId"])
a_species_headers               = np.array(['Sample Name','Assay Name','Description','Protocol REF','Performer','Date','Characteristics [species assay result (VBcv:0000961)]','Term Source Ref','Term Accession Number'])

# @unite-columns:a_species
a_species                       = np.vstack((   col_species_sample_names,\
                                                col_species_assay_names,\
                                                col_species_description,\
                                                col_species_protocolRef,\
                                                col_species_performer,\
                                                col_species_date,\
                                                col_species_result_value,\
                                                col_species_result_termSourceRef,\
                                                col_species_result_accn     ))
#
# rows are flipped to columns
#
a_species = a_species.T

### stack the column headers on top
a_species = np.vstack((a_species_headers,a_species))

# @save:a_species
np.savetxt("../data/isatab/a_species.tsv",      a_species,      delimiter="\t", fmt="%s")


############
# a_IR_WHO #  @@a_IR_WHO    # @done: git commit from home  
############

# Sample Name
# e.g. 'Corrales-El Playón'
col_IR_WHO_sample_names         = np.array([i for i in sample_names])

# Assay Name
# e.g. 'Corrales-El Playón.dr.IR_WHO.lambda-cyhalothrin'
col_IR_WHO_assay_names          = np.array(["PMI.IR_WHO."+i for i in sample_names])
col_IR_WHO_assay_names_final    = []

insecticides = header_to_datacolumn["Insecticide tested"]["raw_dataset_column"]

for i,assay_name in enumerate(col_IR_WHO_assay_names):
    insecticide         = insecticides[i]
    assay_name_final    = assay_name+"."+insecticide
    col_IR_WHO_assay_names_final.append(assay_name_final)

col_IR_WHO_assay_names  = np.array(col_IR_WHO_assay_names_final)

# Protocol REF
# e.g. 'IR_WHO'
col_IR_WHO_protolRef    = np.array(["IR_WHO"]*nrows)

# Performer
# e.g. ''
col_IR_WHO_performer    = np.array([""]*nrows)

# Date
# e.g. ''
col_IR_WHO_date         = np.array([""]*nrows)

# Comment [note]
# e.g. ''
col_IR_WHO_note         = np.array([""]*nrows)

# Parameter Value [group1.insecticidal substance]
# e.g. 'lambda-cyhalothrin'
#col_IR_WHO_insecticide_value = header_to_datacolumn["Insecticide tested"]["raw_dataset_column"]
col_IR_WHO_insecticide_value = header_to_datacolumn["Insecticide tested"]["mapped_ontology_columns"]["ins_term"]

# Term Source Ref
# e.g. 'MIRO'
col_IR_WHO_insecticide_termSourceRef    = np.array(["MIRO"]*nrows)

# Term Accession Number
# e.g. '10000125'
col_IR_WHO_insecticide_accn             = np.array(header_to_datacolumn["Insecticide tested"]["mapped_ontology_columns"]["ins_id"])

# Parameter Value [group1.concentration]
# e.g. '0.05'
col_IR_WHO_concentration_value          = np.array(header_to_datacolumn["Insecticide tested"]["mapped_ontology_columns"]["conc"])

# Unit
# e.g. 'percent'
col_IR_WHO_concentration_unit           = np.array(header_to_datacolumn["Insecticide tested"]["mapped_ontology_columns"]["conc_unit"])

# Term Source Ref
# e.g. 'UO'
col_IR_WHO_concentration_termSourceRef  = np.array(["UO"]*nrows)

# Term Accession Number
# e.g. '0000187'
col_IR_WHO_concentration_accn           = np.array(header_to_datacolumn["Insecticide tested"]["mapped_ontology_columns"]["conc_id"])

# Parameter Value [duration of exposure]
# e.g. '1'
# Unit
# e.g. 'hour'
# Term Source Ref
# e.g. 'UO'
# Term Accession Number
# e.g. '0000032'
col_IR_WHO_durationOfExposure_value = copy(np.array(header_to_datacolumn['Time at which mortality recorded']["raw_dataset_column"]))
col_IR_WHO_durationOfExposure_unit  = copy(col_IR_WHO_durationOfExposure_value)
col_IR_WHO_durationOfExposure_accn  = copy(col_IR_WHO_durationOfExposure_unit)

for i,row in enumerate(col_IR_WHO_durationOfExposure_value):
    if "mins" in row:
        col_IR_WHO_durationOfExposure_value[i]  = row.replace("mins","")
        col_IR_WHO_durationOfExposure_unit[i]   = "minute"
        col_IR_WHO_durationOfExposure_accn[i]   = "0000031"
    elif "hrs" in row:
        col_IR_WHO_durationOfExposure_value[i]  = row.replace("hrs","")
        col_IR_WHO_durationOfExposure_unit[i]   = "hour"
        col_IR_WHO_durationOfExposure_accn[i]   = "0000032"
    else:
        col_IR_WHO_durationOfExposure_value[i]  = ""
        col_IR_WHO_durationOfExposure_unit[i]   = ""
        col_IR_WHO_durationOfExposure_accn[i]   = ""

col_IR_WHO_durationOfExposure_value         = col_IR_WHO_durationOfExposure_value
col_IR_WHO_durationOfExposure_unit          = col_IR_WHO_durationOfExposure_unit
col_IR_WHO_durationOfExposure_termSourceRef = np.array(["UO"]*nrows)
col_IR_WHO_durationOfExposure_accn          = col_IR_WHO_durationOfExposure_accn

# @todo: Additional headers from Unique-values...xlsx: 'Number of replicates for control','Total mosquitoes in all controls','Time at which mortality recorded','Recorded average mortality in treatments (%)','Recorded average mortality in controls (%)'

# Raw Data File
# e.g. 'p_IR_WHO.txt'
col_IR_WHO_rawDataFile = np.array(['p_IR_WHO.txt']*nrows)

##  <<<@DONE: do header numbers>>>
## <<<<<<<@todo: add new headers in "Comment []" form, after browsing: Unique-values...xlsx
a_IR_WHO_headers = np.array(['Sample Name','Assay Name','Protocol REF','Performer','Date','Comment [note]','Parameter Value [group1.insecticidal substance]','Term Source Ref','Term Accession Number','Parameter Value [group1.concentration]','Unit','Term Source Ref','Term Accession Number','Parameter Value [duration of exposure]','Unit','Term Source Ref','Term Accession Number','Raw Data File'])

a_IR_WHO = np.array([   col_IR_WHO_sample_names,\
                        col_IR_WHO_assay_names_final,\
                        col_IR_WHO_protolRef,\
                        col_IR_WHO_performer,\
                        col_IR_WHO_date,\
                        col_IR_WHO_note,\
                        col_IR_WHO_insecticide_value,\
                        col_IR_WHO_insecticide_termSourceRef,\
                        col_IR_WHO_insecticide_accn,\
                        col_IR_WHO_concentration_value,\
                        col_IR_WHO_concentration_unit,\
                        col_IR_WHO_concentration_termSourceRef,\
                        col_IR_WHO_concentration_accn,\
                        col_IR_WHO_durationOfExposure_value,\
                        col_IR_WHO_durationOfExposure_unit,\
                        col_IR_WHO_durationOfExposure_termSourceRef,\
                        col_IR_WHO_durationOfExposure_accn,\
                        col_IR_WHO_rawDataFile  ])

# rows are flipped to columns
a_IR_WHO = a_IR_WHO.T

# filter out non IR_BA rows 
a_IR_WHO_filtered = filter_assay_rows( a_IR_WHO, "WHO test kit_adults" )  # @DONE: check length of filtered WHO assay is 

# stack the column headers on top
a_IR_WHO = np.vstack(( a_IR_WHO_headers, a_IR_WHO_filtered ))

# @save:a_IR_WHO
np.savetxt("../data/isatab/a_IR_WHO.tsv",      a_IR_WHO,      delimiter="\t", fmt="%s")

############
# p_IR_WHO #          @@p_IR_WHO
############
#   <<<<< @todo: need to filter to include only WHO_testkit rows >>>>>
#   <<<<< @todo: sanity check the code up till p_IR_BA, as they were done ad hoc from home... >>>>>

# Assay Name
# e.g.: Corrales-El Playón.dr.IR_WHO.deltamethrin
col_p_IR_WHO_assay_names = np.array(["PMI.IR_WHO."+i for i in sample_names])
col_p_IR_WHO_assay_names = col_IR_WHO_assay_names

# Phenotype Name
# e.g.: Mortality percentage:100, 0.05% deltamethrin
col_p_IR_WHO_mortalityPercentage    = []
mortalityPercentages                = header_to_datacolumn['Calculated  average mortality adjusted for control  (%)']["raw_dataset_column"]  # e.g. NR, 0.987, 42.58%
mortalityPercentages_no_unit        = []

for i,mortality in enumerate(mortalityPercentages):
    # e.g. Mortality percentage:100, 0.05% deltamethrin
    mortality_noUnit    = mortality.replace("%","")
    insecticide         = col_IR_WHO_insecticide_value[i]
    concentration_value = col_IR_WHO_concentration_value[i]
    concentration_unit  = col_IR_WHO_concentration_unit[i]

    if concentration_unit == "percent":
        concentration_unit = "%"

    col_p_IR_WHO_mortalityPercentage.append("mortality percentage:"+mortality_noUnit+", "+concentration_value+" "+concentration_unit+" "+insecticide)  # e.g. Mortality percentage:100, 0.05% lambda-cyhalothrin

    mortalityPercentages_no_unit.append(mortality_noUnit)

# @DONE: see if the Assay name and Phenotype name columns look as they should vs. fonseca

# Observable
# e.g.: insecticide resistance
col_p_IR_WHO_observable_value         = np.array(["insecticide resistance"]*nrows)

# Term Source Ref
# e.g.: MIRO
col_p_IR_WHO_observable_termSourceRef = np.array(["MIRO"]*nrows)

# Term Accession Number
# e.g.: 00000021
col_p_IR_WHO_observable_accn          = np.array(["00000021"]*nrows)

# Attribute
# e.g.: mortality rate
col_p_IR_WHO_attribute_value          = np.array(["mortality rate"]*nrows)

# Term Source Ref
# e.g.: VBcv
col_p_IR_WHO_attribute_termSourceRef  = np.array(["VBcv"]*nrows)

# Term Accession Number
# e.g.: 0000703
col_p_IR_WHO_attribute_accn           = np.array(["0000703"]*nrows)

# Value
# e.g.: 100
col_p_IR_WHO_value_value              = mortalityPercentages_no_unit

# Term Source Ref
# e.g.: ""
col_p_IR_WHO_value_termSourceRef      = np.array([""]*nrows)

# Term Accession Number
# e.g.: ""
col_p_IR_WHO_value_accn               = np.array([""]*nrows)

# Unit
# e.g.: percent
col_p_IR_WHO_value_unit               = np.array(["percent"]*nrows)

# Term Source Ref
# e.g.: UO
col_p_IR_WHO_value_unit_termSourceRef = np.array(["UO"]*nrows)

# Term Accession Number
# e.g.: 0000187
col_p_IR_WHO_value_unit_accn          = np.array(["0000187"]*nrows)

p_IR_WHO_headers = np.array([   'Assay Name',\
                                'Phenotype Name',\
                                'Observable',\
                                'Term Source Ref',\
                                'Term Accession Number',\
                                'Attribute',\
                                'Term Source Ref',\
                                'Term Accession Number',\
                                'Value',\
                                'Term Source Ref',\
                                'Term Accession Number',\
                                'Unit',\
                                'Term Source Ref',\
                                'Term Accession Number'     ]) 

p_IR_WHO = np.array([           col_p_IR_WHO_assay_names,
                                col_p_IR_WHO_mortalityPercentage,
                                col_p_IR_WHO_observable_value,
                                col_p_IR_WHO_observable_termSourceRef,
                                col_p_IR_WHO_observable_accn,
                                col_p_IR_WHO_attribute_value,
                                col_p_IR_WHO_attribute_termSourceRef,
                                col_p_IR_WHO_attribute_accn,
                                col_p_IR_WHO_value_value,
                                col_p_IR_WHO_value_termSourceRef,
                                col_p_IR_WHO_value_accn,
                                col_p_IR_WHO_value_unit,
                                col_p_IR_WHO_value_unit_termSourceRef,
                                col_p_IR_WHO_value_unit_accn    ])

# rows are flipped to columns
p_IR_WHO             = p_IR_WHO.T


#pdb.set_trace()  # @0950

# filter out non IR_BA rows 
p_IR_WHO_filtered    = filter_assay_rows( p_IR_WHO, "WHO test kit_adults" )  # @DONE: check length of filtered WHO assay is 

# stack the column headers on top
p_IR_WHO             = np.vstack(( p_IR_WHO_headers, p_IR_WHO_filtered ))

# @save:p_IR_WHO
np.savetxt("../data/isatab/p_IR_WHO.tsv",      p_IR_WHO,      delimiter="\t", fmt="%s")


############
# a_IR_BA  #    @@a_IR_BA
############

#   <<<<< @todo: need to filter to include only WHO_testkit rows >>>>>
#   <<<<< @todo: sanity check the code up till p_IR_BA, as they were done ad hoc from home... >>>>>
col_IR_BA_sample_names = col_IR_WHO_sample_names
col_IR_BA_assay_names  = col_IR_WHO_assay_names

col_IR_BA_assay_names_final                 = [i.replace("IR_WHO","IR_BA") for i in col_IR_BA_assay_names]
col_IR_BA_protolRef                         = copy([i.replace("IR_WHO","IR_BA") for i in col_IR_BA_protolRef])
col_IR_BA_performer                         = copy(col_IR_WHO_performer)
col_IR_BA_date                              = copy(col_IR_WHO_date)
col_IR_BA_note                              = copy(col_IR_WHO_note)
col_IR_BA_insecticide_value                 = copy(col_IR_WHO_insecticide_value)
col_IR_BA_insecticide_termSourceRef         = copy(col_IR_WHO_insecticide_termSourceRef)
col_IR_BA_insecticide_accn                  = copy(col_IR_WHO_insecticide_accn)
col_IR_BA_concentration_value               = copy(col_IR_WHO_concentration_value)
col_IR_BA_concentration_unit                = copy(col_IR_WHO_concentration_unit)
col_IR_BA_concentration_termSourceRef       = copy(col_IR_WHO_concentration_termSourceRef)
col_IR_BA_concentration_accn                = copy(col_IR_WHO_concentration_accn)
col_IR_BA_durationOfExposure_value          = copy(col_IR_WHO_durationOfExposure_value)
col_IR_BA_durationOfExposure_unit           = copy(col_IR_WHO_durationOfExposure_unit)
col_IR_BA_durationOfExposure_accn           = copy(col_IR_WHO_durationOfExposure_accn)
col_IR_BA_durationOfExposure_termSourceRef  = copy(col_IR_WHO_durationOfExposure_termSourceRef)
col_IR_BA_rawDataFile                       = copy(col_IR_WHO_rawDataFile)

a_IR_BA_headers = np.array([    'Sample Name',\
                                'Assay Name',\
                                'Protocol REF',\
                                'Performer',\
                                'Date',\
                                'Comment [note]',\
                                'Parameter Value [group1.insecticidal substance]',\
                                'Term Source Ref',\
                                'Term Accession Number',\
                                'Parameter Value [group1.concentration]',\
                                'Unit',\
                                'Term Source Ref',\
                                'Term Accession Number',\
                                'Parameter Value [duration of exposure]',\
                                'Unit',\
                                'Term Source Ref',\
                                'Term Accession Number',\
                                'Raw Data File' ])

a_IR_BA = np.array([    col_IR_BA_sample_names,\
                        col_IR_BA_assay_names_final,\
                        col_IR_BA_protolRef,\
                        col_IR_BA_performer,\
                        col_IR_BA_date,\
                        col_IR_BA_note,\
                        col_IR_BA_insecticide_value,\
                        col_IR_BA_insecticide_termSourceRef,\
                        col_IR_BA_insecticide_accn,\
                        col_IR_BA_concentration_value,\
                        col_IR_BA_concentration_unit,\
                        col_IR_BA_concentration_termSourceRef,\
                        col_IR_BA_concentration_accn,\
                        col_IR_BA_durationOfExposure_value,\
                        col_IR_BA_durationOfExposure_unit,\
                        col_IR_BA_durationOfExposure_termSourceRef,\
                        col_IR_BA_durationOfExposure_accn,\
                        col_IR_BA_rawDataFile ])

# rows are flipped to columns
a_IR_BA             = a_IR_BA.T

# filter out non IR_BA rows 
a_IR_BA_filtered    = filter_assay_rows( a_IR_BA, "CDC bottle_adults" )  # @DONE: check length of filtered WHO assay is 

# stack the column headers on top
a_IR_BA             = np.vstack((a_IR_BA_headers,a_IR_BA_filtered))

# @save:a_IR_BA
np.savetxt("../data/isatab/a_IR_BA.tsv",      a_IR_BA,      delimiter="\t", fmt="%s")

############
# p_IR_BA  #   @done: git commit from home
############

# @todo: put these rows of data (col_p_IR_*) stacked then headers then filter using the filter_assay_rows() functions                @latest

col_p_IR_BA_assay_names              = copy(col_p_IR_WHO_assay_names)
col_p_IR_BA_mortalityPercentage      = copy(col_p_IR_WHO_mortalityPercentage)
col_p_IR_BA_observable_value         = copy(col_p_IR_WHO_observable_value)
col_p_IR_BA_observable_termSourceRef = copy(col_p_IR_WHO_observable_termSourceRef)
col_p_IR_BA_observable_accn          = copy(col_p_IR_WHO_observable_accn)
col_p_IR_BA_attribute_value          = copy(col_p_IR_WHO_attribute_value)
col_p_IR_BA_attribute_termSourceRef  = copy(col_p_IR_WHO_attribute_termSourceRef)
col_p_IR_BA_attribute_accn           = copy(col_p_IR_WHO_attribute_accn)
col_p_IR_BA_value_value              = copy(col_p_IR_WHO_value_value)
col_p_IR_BA_value_termSourceRef      = copy(col_p_IR_WHO_value_termSourceRef)
col_p_IR_BA_value_accn               = copy(col_p_IR_WHO_value_accn)
col_p_IR_BA_value_unit               = copy(col_p_IR_WHO_value_unit)
col_p_IR_BA_value_unit_termSourceRef = copy(col_p_IR_WHO_value_unit_termSourceRef)
col_p_IR_BA_value_unit_accn          = copy(col_p_IR_WHO_value_unit_accn)

# make compatible w/ p_IR_BA_* format 
p_IR_BA_headers = np.array([   'Assay Name',\
                                'Phenotype Name',\
                                'Observable',\
                                'Term Source Ref',\
                                'Term Accession Number',\
                                'Attribute',\
                                'Term Source Ref',\
                                'Term Accession Number',\
                                'Value',\
                                'Term Source Ref',\
                                'Term Accession Number',\
                                'Unit',\
                                'Term Source Ref',\
                                'Term Accession Number'     ]) 

p_IR_BA = np.array([    col_p_IR_BA_assay_names,\
                        col_p_IR_BA_mortalityPercentage,\
                        col_p_IR_BA_observable_value,\
                        col_p_IR_BA_observable_termSourceRef,\
                        col_p_IR_BA_observable_accn,\
                        col_p_IR_BA_attribute_value,\
                        col_p_IR_BA_attribute_termSourceRef,\
                        col_p_IR_BA_attribute_accn,\
                        col_p_IR_BA_value_value,\
                        col_p_IR_BA_value_termSourceRef,\
                        col_p_IR_BA_value_accn,\
                        col_p_IR_BA_value_unit,\
                        col_p_IR_BA_value_unit_termSourceRef,\
                        col_p_IR_BA_value_unit_accn ])

# rows are flipped to columns
p_IR_BA             = p_IR_BA.T

# filter out non p_IR_BA rows 
p_IR_BA_filtered    = filter_assay_rows( p_IR_BA, "CDC bottle_adults" )  # @DONE: check length of filtered WHO assay is 

# stack the column headers on top
p_IR_BA             = np.vstack(( p_IR_BA_headers, p_IR_BA_filtered ))

# @save:a_IR_BA
np.savetxt("../data/isatab/p_IR_BA.tsv",      p_IR_BA,      delimiter="\t", fmt="%s")


############################################################################

################
################
# save to file #
################
################

# #
# # check if necessary paths exist, else make them
# #
# # @DONE: check if file exists else create it
# check_path_exists_else_make_it("../data/")
# check_path_exists_else_make_it("../data/ontologies")
# check_path_exists_else_make_it("../data/raw")
# check_path_exists_else_make_it("../data/isatab")

# #
# # write ISA-Tab files  
# #
# np.savetxt("../data/isatab/s_samples.tsv",      s_samples,      delimiter="\t", fmt="%s")
# np.savetxt("../data/isatab/a_collection.tsv",   a_collection,   delimiter="\t", fmt="%s")
# np.savetxt("../data/isatab/a_species.tsv",      a_species,      delimiter="\t", fmt="%s")

# #
# # @DONE: check that columns of these two sheets match against Fonseca format
# #
# np.savetxt("../data/isatab/a_IR_WHO.tsv",       a_IR_WHO,       delimiter="\t", fmt="%s")
# np.savetxt("../data/isatab/a_IR_BA.tsv",        a_IR_BA,        delimiter="\t", fmt="%s")

# #
# # @DONE: check that columns of these two sheets match against Fonseca format
# #
# np.savetxt("../data/isatab/p_IR_WHO.tsv",       p_IR_WHO,        delimiter="\t", fmt="%s")
# np.savetxt("../data/isatab/p_IR_BA.tsv",        p_IR_BA,        delimiter="\t", fmt="%s")

print("--- %s seconds ---" % (time.time() - start_time))