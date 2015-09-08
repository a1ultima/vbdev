
""" Generates foreign keys mapping: Raw PMI values -to- VectorBase-friendly ontology terms, in the format of a "python dict", which are then stored as a pickle (.p) files

For example /home/ab108/0VB/isatab/presidentsMalariaInitiative/data/ontologies/country_to_gaz.txt looks like:

    Country Characteristics [Collection site (VBcv:0000831)]        Term Source Ref Term Accession Number
    Ghana   Ghana   GAZ     00000908
    Liberia Liberia GAZ     00000911
    Madagascar      Madagascar      GAZ     00001108
    Angola  Angola  GAZ     00001095
    Ethiopia        Ethiopia        GAZ     00000567
    Zambia  Zambia  GAZ     00001107
    Zimbabwe        Zimbabwe        GAZ     00001106
    Malawi  Malawi  GAZ     00001105
    Democratic Republic of the Congo        Democratic Republic of the Congo        GAZ     00001086
    Senegal Senegal GAZ     00000913
    Rwanda  Rwanda  GAZ     00001087
    Mozambique      Mozambique      GAZ     00001100
    Uganda  Uganda  GAZ     00001102
    Mali    Mali    GAZ     00000584
    Burundi Burundi GAZ     00001090
    Kenya   Kenya   GAZ     00001101
    Nigeria Nigeria GAZ     00000912
    Tanzania, United Republic of    Tanzania        GAZ     00001103
    Burkina Faso    Burkina Faso    GAZ     00000905
    Benin   Benin   GAZ     00000904

And there are 8 of these (in /home/ab108/0VB/isatab/presidentsMalariaInitiative/data/ontologies):

    country_to_gaz.txt
    deleteMeAfterThis.txt
    endMonth_to_numberMonth.txt
    insecticide_to_ontoTerms.txt
    speciesTested_to_miro.txt
    stageTested_to_ontoTermsAdult.txt
    startMonth_to_numberMonth.txt
    yearStartEnd_to_numStartEnd.txt

Which would then be stored as "pickled python dicts" (.p) files:

    country_to_gazTerm_gazId.p
    insecticide_to_ontoTerms.p
    species_to_miroTerm_miroId.p
    yearStartToEnd_yearStart_yearEnd.p
    startMonth_to_numMonth.p
    stageTested_ontoTermAdult_termSourceRef_accnNum.p
    endMonth_to_numMonth.p


"""

#############
# Imports   #
#############

import csv 
import numpy as np
import pdb
import pickle

#############
# Functions #
#############

def pickle_dump(objName,fName):
    f = open('../data/ontologies/'+fName,'wb')
    pickle.dump(objName,f)
    f.close()

#######
# Run #
#######

print "Reading in foreign keys mapping from: Raw PMI values -to- VectorBase-friendly ontology terms..."

#
# country to: gazTerm, gazId
#
reader       = csv.reader(open("../data/ontologies/country_to_gaz.txt","rb"),delimiter="\t")
country_rows = list(reader)
country_to_gazTerm_gazId = {}

for row in country_rows:
    # e.g. row: Ghana   Ghana   GAZ 00000908
    country = row[0]
    gazTerm = row[1]
    source  = row[2]  # always: GAZ
    gazId   = row[3]
    country_to_gazTerm_gazId[country]={'gazTerm':gazTerm,'source':source,'gazId':gazId}

#
# speciesTested to: miroTerm, miroId
#
reader      = csv.reader(open("../data/ontologies/speciesTested_to_miro.txt","rb"),delimiter="\t")
print "\t/home/ab108/0VB/isatab/presidentsMalariaInitiative/data/ontologies/data/ontologies/speciesTested_to_miro.txt"
species_rows= list(reader)
species_to_miroTerm_miroId = {}

for row in species_rows:
    # e.g. An. pharoensis  Anopheles pharoensis    MIRO    40003564
    species  = row[0]
    miroTerm = row[1]
    source   = row[2]
    miroId   = row[3]
    species_to_miroTerm_miroId[species]={'miroTerm':miroTerm,'miroId':miroId}

#
# insecticide tested to: miroTerms, concentrations, etc. 
#
reader  = csv.reader(open("../data/ontologies/insecticide_to_ontoTerms.txt","rb"),delimiter="\t")
print "\t/home/ab108/0VB/isatab/presidentsMalariaInitiative/data/ontologies/data/ontologies/insecticide_to_ontoTerms.txt"
insecticide_rows = list(reader)  # e.g. Malathion 5.0%  malathion   MIRO    10000065    5   percent UO  0000187
insecticide_to_ontoTerms = {}  # e.g. country_to_gazTerm_gazId
for row in insecticide_rows:
    # e.g. Malathion 5.0%   malathion   MIRO    10000065    5   percent UO  0000187
    ins_tested    = row[0]  # insecticide_tested
    ins_term      = row[1]  # miroTerm = row[1]
    ins_source    = row[2]  # miroId   = row[2]
    ins_id        = row[3]  # miroId   = row[2]
    conc          = row[4]  # e.g. 0.05
    conc_unit     = row[5]  # e.g. aeg / %
    conc_source   = row[6]  # always UO
    conc_id       = row[7]  # e.g. 0000187
    insecticide_to_ontoTerms[ins_tested]={'ins_term':ins_term,'ins_source':ins_source,'ins_id':ins_id,'conc':conc,'conc_unit':conc_unit,'conc_source':conc_source,'conc_id':conc_id}  # insecticide_to_ontoTerms[insecticide]={'miroTerm':miroTerm,'miroId':miroId}

#
# startMonth_to_numberMonth
#
reader  = csv.reader(open("../data/ontologies/startMonth_to_numberMonth.txt","rb"),delimiter="\t")
print "\t/home/ab108/0VB/isatab/presidentsMalariaInitiative/data/ontologies/data/ontologies/startMonth_to_numberMonth.txt"
startMonth_rows = list(reader)  # e.g. species
startMonth_to_numMonth = {}  # e.g. country_to_gazTerm_gazId

for row in startMonth_rows:
    # e.g. 
    startMonth = row[0]  # species
    numMonth = row[1]  # miroTerm = row[1]
    startMonth_to_numMonth[startMonth] = numMonth  # startMonth_to_numMonth[startMonth]={'miroTerm':miroTerm,'miroId':miroId}

#
# endMonth_to_numberMonth
#
reader  = csv.reader(open("../data/ontologies/endMonth_to_numberMonth.txt","rb"),delimiter="\t")
print "\t/home/ab108/0VB/isatab/presidentsMalariaInitiative/data/ontologies/data/ontologies/endMonth_to_numberMonth.txt"
endMonth_rows = list(reader)  # e.g. species
endMonth_to_numMonth = {}  # e.g. country_to_gazTerm_gazId

for row in endMonth_rows:
    endMonth = row[0]  # species
    numMonth = row[1]  # miroTerm = row[1]
    endMonth_to_numMonth[endMonth]= numMonth

#
# year range to yearStart, yearEnd
#
reader  = csv.reader(open("../data/ontologies/yearStartEnd_to_numStartEnd.txt","rb"),delimiter="\t")
print "\t/home/ab108/0VB/isatab/presidentsMalariaInitiative/data/ontologies/data/ontologies/yearStartEnd_to_numStartEnd.txt"
yearStartToEnd_rows = list(reader)  # e.g. species
yearStartToEnd_yearStart_yearEnd = {}  # e.g. country_to_gazTerm_gazId
for row in yearStartToEnd_rows:
    yearStartToEnd = row[0]  # species
    yearStart = row[1]  # miroTerm = row[1]
    yearEnd = row[2]  # miroId   = row[2]
    yearStartToEnd_yearStart_yearEnd[yearStartToEnd]={'yearEnd':yearEnd,'yearStart':yearStart}  # yearStartToEnd_yearStart_yearEnd[yearStartToEnd]={'miroTerm':miroTerm,'miroId':miroId}

#
# Stage tested 
#
reader  = csv.reader(open("../data/ontologies/stageTested_to_ontoTermsAdult.txt","rb"),delimiter="\t")
print "\t/home/ab108/0VB/isatab/presidentsMalariaInitiative/data/ontologies/data/ontologies/stageTested_to_ontoTermsAdult.txt"
stageTested_rows = list(reader)  # e.g. species
stageTested_ontoTermAdult_termSourceRef_accnNum = {}  # e.g. country_to_gazTerm_gazId
for row in stageTested_rows:
    stageTested = row[0]  # species
    ontoTermAdult = row[1]  # miroTerm = row[1]
    termSourceRef = row[2]  # miroId   = row[2]
    accn = row[3]
    stageTested_ontoTermAdult_termSourceRef_accnNum[stageTested]={'ontoTermAdult':ontoTermAdult,'termSourceRef':termSourceRef,'accn':accn}  

#
# Sanity Checking
#
# print "\tNumbers of values mapped in each foreign key..."
# print "\t\tspecies_to_miroTerm_miroId "+str(len(species_to_miroTerm_miroId.keys()))
# print "\t\tcountry_to_gazTerm_gazId "+str(len(country_to_gazTerm_gazId.keys()))
# print "\t\tinsecticide_to_ontoTerms "+str(len(insecticide_to_ontoTerms.keys()))
# print "\t\tspecies_to_miroTerm_miroId"+str(len(species_to_miroTerm_miroId))
# print "\t\tyearStartToEnd_yearStart_yearEnd "+str(len(yearStartToEnd_yearStart_yearEnd.keys()))
# print "\t\tstartMonth_to_numMonth "+str(len(startMonth_to_numMonth.keys()))
# print "\t\tendMonth_to_numMonth "+str(len(endMonth_to_numMonth.keys()))
# print "\t\tstageTested_to_ontoTermsAdult"+str(len(stageTested_ontoTermAdult_termSourceRef_accnNum.keys()))

#
# Dump the ontology keys as dictionaries to be loaded by the ISA-Tab generator
#
pickle_dump(species_to_miroTerm_miroId,"species_to_miroTerm_miroId.p")                  # 1
pickle_dump(country_to_gazTerm_gazId,"country_to_gazTerm_gazId.p")                      # 2
pickle_dump(insecticide_to_ontoTerms,"insecticide_to_ontoTerms.p")                      # 3
pickle_dump(species_to_miroTerm_miroId,"species_to_miroTerm_miroId.p")                  # 4
pickle_dump(yearStartToEnd_yearStart_yearEnd,"yearStartToEnd_yearStart_yearEnd.p")      # 5
pickle_dump(startMonth_to_numMonth,"startMonth_to_numMonth.p")                          # 6
pickle_dump(endMonth_to_numMonth,"endMonth_to_numMonth.p")                              # 7
pickle_dump(stageTested_ontoTermAdult_termSourceRef_accnNum,"stageTested_ontoTermAdult_termSourceRef_accnNum.p")  # 8
