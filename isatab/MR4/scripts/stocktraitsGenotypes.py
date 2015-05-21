
""" README: ./scripts/stocktraitsGenotypes.py


Description:

    Converts raw data of inversion karyotypes for 2L and 2R available in "Stock traits for Vectorbase continued.xlsx" from https://www.ebi.ac.uk/panda/jira/browse/VB-2088 into ISA-Tab friendly format. The relevant region of the .xlsx is: columns D:F, rows 2:23. These are copied from the .xlsx formatted into python list, then pasted into this script: CTRL+F: "@paste"

    Example of desired format for entry in g_karyotypes.txt ISA-Tab format is found in the Neafsey ISA-Tab: https://docs.google.com/spreadsheet/ccc?key=0AiHJk7ao3qFUdEpWM2NYWWhZa0tmaDlyWll3MUpITGc&usp=sheets_web#gid=8    

Usage:

    From shell: 

        cd /home/ab108/0VB/isatab/MR4/scripts/
        python ./stocktraitsGenotypes.py > ../data/g_karyotypes.txt

"""



#########################################
# DATA COPY-PASTED FROM STOCKTRAITS     #  @paste
#########################################

# column 'C', row '2:23'
names  =['MRA-121','MRA-913','MRA-114','MRA-860','MRA-763','MRA-765','MRA-186','MRA-112','MRA-115','MRA-762','MRA-105','MRA-861','MRA-111','MRA-334','MRA-698','MRA-594','MRA-339','MRA-764','MRA-856','MRA-1155','MRA-803','MRA-1156']

# column 'D:F', row '1'
headers=['2La (PCR)','2Rb (PCR)','2Ro (PCR)']
headers=[header.replace(' (PCR)','') for header in headers] # reformat

# column 'D:F', row '2:23'
columns=[['+/La','+/b','+/+'],['+/La','+/+','+/+'],['+/+','+/b','+/+'],['La/La','b/b','+/+'],['La/La','+/+','+/+'],['La/La','--','+/+'],['+/La','+/+','+/+'],['+/La','+/b','+/+'],['+/La','b/b','+/+'],['+/+','+/b','+/+'],['+/+','+/+','+/+'],['+/La','b/b','+/+'],['+/La','+/+','+/+'],['+/La','+/+','+/+'],['+/+','+/+','+/+'],['+/La','+/+','+/+'],['La/La','+/+','+/+'],['La/La','+/+','+/+'],['La/La','+/+','+/+'],['+/+','',''],['La/La','','o/o'],['La/La','','o/o']]



############################################################################################################
# Convert StockTraits.xls format of inversion genotypes into ISA-Tab friendly g_karyotypes.txt format
############################################################################################################

L2_vec = []
R2_vec = []

for i, (L2a, R2b, R2o) in enumerate(columns):

    #print i, L2a, R2b, R2o

    # Reformat the 2L chromosome inversion karyotypes
    if L2a == '+/La':
        L2_vec.append('2La/+')
    elif L2a == '+/+':
        L2_vec.append('2L+/+')
    elif L2a == 'La/La':
        L2_vec.append('2La/a')
    # else:
    #     #raise Exception('Unexpected karyotype!: '+str(L2a)+' for: '+str(names[i]))
    #     print 'Unexpected karyotype!: '+str(L2a)+' for: '+str(names[i])
    #     pass


    # Reformat the 2R chromosome inversion karyotypes
    R2_slot1 = ''
    R2_slot2 = ''

    if R2b == '+/b':
        R2_slot1 += 'b'
    elif R2b == '+/+':
        pass
    elif R2b == 'b/b':
        R2_slot1 += 'b'
        R2_slot2 += 'b'
    # else:
    #     #raise Exception('Unexpected karyotype!: '+str(R2b)+' for: '+str(names[i]))
    #     print 'Unexpected karyotype!: '+str(R2b)+' for: '+str(names[i])
    #     pass

    if R2o == '+/o':
        R2_slot1 += 'o'
    elif R2o == '+/+':
        pass
    elif R2o == 'o/o':
        R2_slot1 += 'o'
        R2_slot2 += 'o'
    # else:
        # #raise Exception('Unexpected karyotype!: '+str(R2o)+' for: '+str(names[i]))
        # print 'Unexpected karyotype!: '+str(R2o)+' for: '+str(names[i])
        # pass


    if R2_slot1 == '':
        R2_slot1 = '+'
    if R2_slot2 == '':
        R2_slot2 = '+'

    R2 = '2R'+R2_slot1+'/'+R2_slot2

    R2_vec.append(R2)


###########################################
#  Generate Columns of ISA-Tab ready data #
###########################################

headers = 'Assay Name'+'\t'+'Genotype Name'+'\t'+'Description'+'\t'+'Type'+'\t'+'Term Source REF'+'\t'+'Term Accession Number'+'\t'+'Characteristics [inversion (SO:1000036)]'+'\t'+'Term Source REF'+'\t'+'Term Accession Number'+'\t'+'Characteristics [chromosome_arm (SO:0000105)]'+'\t'+'Term Source REF'+'\t'+'Term Accession Number'

print headers

assay_names = []

for i,name in enumerate(names):

    row1 = name+'.ka1'+'\t'+L2_vec[i]+'\t'+'inversion: '+L2_vec[i]+'\t'+'paracentric_inversion'+'\t'+'SO'+'\t'+'1000047'+'\t'+L2_vec[i]+'\t'+''+'\t'+''+'\t'+'2L'+'\t'+''+'\t'+''
    row2 = name+'.ka1'+'\t'+R2_vec[i]+'\t'+'inversion: '+R2_vec[i]+'\t'+'paracentric_inversion'+'\t'+'SO'+'\t'+'1000047'+'\t'+R2_vec[i]+'\t'+''+'\t'+''+'\t'+'2R'+'\t'+''+'\t'+''
    
    assert headers.count('\t') == row1.count('\t'), 'Headers and row1 have different numbers of columns!'
    assert headers.count('\t') == row2.count('\t'), 'Headers and row2 have different numbers of columns!'

    print row1
    print row2

###########
# Testing #
###########

    #----------------------------------------------------------------------------------------------------------------
    # TEST: MRA names from "stocktraits continued.xlsx" are present in the ISA-Tab s_samples.txt names
    #----------------------------------------------------------------------------------------------------------------

sample_names = ['MRA-105','MRA-111','MRA-112','MRA-114','MRA-115','MRA-121','MRA-126','MRA-128','MRA-130','MRA-139','MRA-186','MRA-314','MRA-334','MRA-339','MRA-489','MRA-493','MRA-594','MRA-698','MRA-700','MRA-726','MRA-728','MRA-729','MRA-730-CDC','MRA-734','MRA-735','MRA-762','MRA-763','MRA-764','MRA-765','MRA-803','MRA-804','MRA-854','MRA-856','MRA-860','MRA-861','MRA-862','MRA-863','MRA-864','MRA-891','MRA-893','MRA-913','MRA-978','MRA-1027','MRA-1154','MRA-1155','MRA-1156','MRA-1163','MRA-1164','MRA-1165','MRA-1166','MRA-1167']

for name in names:

    assert name in sample_names, name +'from "Stock traits for Vectorbase continued.xlsx" not present in s_samples.txt' 

