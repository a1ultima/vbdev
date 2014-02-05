from cogent.db.ensembl import HostAccount, Genome, Species

# anopheles_albimanus_core_1312_73_1

# example human antisense gene: 'ENST00000430895'

account = None
genome = Genome(Species='Human',Release='73',account=account)

#______________________________________________________
# TEST 1: absolute coords or relative to -/+ strand?
#______________________________________________________

# # A known antisense
# gene1 = genome.getGeneByStableId('ENSG00000234484')       # Coordinate(Human,chro...,6,133073813-133075090,1)

# # A sense gense upstream from ^ but on - strand
# gene2 = genome.getGeneByStableId('ENSG00000146399')       # Coordinate(Human,chro...,6,133073813-133075090,1)

# # A sense gene downstream from ^ but on + strand
# gene3 = genome.getGeneByStableId('ENSG00000271488')

# # Andy's predictions: if absolute coords are used, then gene2 < gene1 < gene3 == True
# print gene2.Location < gene1.Location < gene3.Location # Indeed, it seems absolute coords are used... YAY! 

Species.amendSpecies('anopheles albimanus','anopheles_albimanus')
account = HostAccount('localhost','vbuser','Savvas')
genome = Genome(Species='anopheles albimanus',Release='73',account=account)
gene = genome.getGeneByStableId('AALB004251')

#gene = genome.getGenesMatching('a')

#______________________________________________________
# TEST 2: Check if 2kb really is 2kb and not +/- 1 bp
#______________________________________________________

# Here's the printout of part of the sample_data dict loaded from pickle
# On sublime: Double click on "untruncated" then look at how many chars are selected
"""
{'truncated': '', 'untruncated': ''},
 {'truncated': 'TAGGCAATCGGTGTCCCATTTTTTGTTGCAGCAACAGCAGCAGCAGCAGCCTCAGCAGCCGCAACAGCAATCATTCCAGCAGCAGCAACAAGCACCCCCAGATCTGCAACAGAGCCTGCAGCATTCGCGGCATCAACTGCAGCTACAGCAAAACCAGTTGCAGCAACAGTTACAGCAGCAACAACAACAACATCGTACGCAGCATCAGCAATTGCAACAATTGCAGGCTCACCAGCAACATCATCCTGCATCACAGCAACAACCTCAGCAGACTCAGCTCATGGTACAAACGATGTCGGAGCCGGATTTGTTTGGTCCTAAGATAGCG',
  'untruncated': 'ATCGACTGTGCTGAAGTCGACGATTGCGTTCCTGAAGAGCCACAACGAGATCGCGGTACGGTCACGGGTACACGAGATCCAGACGGATTGGAAGCCATCGTTCCTGTCGAACGAAGAGTTTACGCATCTGATTCTGGAAGCACTGGACGGGTTCATCATCGTGTTCTCTTCCACCGGCAGAGTATTTTATGCGTCCGAAAGTATTACCTCCTTGTTGGGCCACTTACCGGTTAGTAGTGATTGGATGGTGTGACGTAGTATTTGCTAACGGTTTTGATGTTTTTTGCTGTTGTGATTCCAGAGTGATCTGCTCAACATGACCGTGTACGATATGGTGTACGAGGACGATCAGAATGATCTGTATCACATACTGCTTAACCCGACCACGATCGTGGATCCGCTTCAGACGGGCATTAGCCGAGAGAATCAGGTAACCTTTTCCTGTTACATCAAACGTGGAACGGTTGACTATCGGGCGGATGTGTCGTACGAGCATGTCCAGTTCACTGGGTACTTCCGTGAGTATACTAAGGCATGAAGAAGAAGGTCATGAAAGAGGTAAACATGGTGCTGAATGATGGCTCTTCGTAGGTAGTGACGTCGATACGGAATCTCTTATGACAACATCTCGGTTCAGCGGCTATACGAGCGATGCCGACTCTCGACTCGTGTTTGTAGGTACGGGTCGGCTGCAGACCCCGCAACTGATCCGCGAGATGTCGATCGTGGATAACACGAAAAGTGAGTTCACCTCCCGGCACAGCCTCGAATGGAAGTTCCTGTTTCTCGATCACCGGGCACCACCTATCATCGGCTACCTCCCGTTTGAAGTGCTGGGAACTTCCGGGTATGACTATTATCATTTTGACGATCTGGAAAAGGTGGTGGCTTGTCATGAAGCATGTAAGCTCTCACAGTATTACTACTACTGCCGTGCCCTCGTTCGCTACGGCGATCGCCTCTCCAAACGCATGGCCTGGCTTCACCGTACTCTTTCTCTTTTTCTTTACCTTCCCGTTTGTGTCTCTCTACTACTCTCGGCGTTGCTCGCCGTTTTGACGCTTATAGTGATGCAGAAAGGTGAGGGCACTTCGTGCTACTATCGATTCCTGACTAAAGGGCAGCAGTGGATCTGGCTGCAGACACGGTTCTACATCACCTACCATCAGTGGAACTCGAAGCCGGAGTTCGTGGTGTGTACGCATCGAGTGGTAAGCTATGCGGATGTGATGAAACAGATGCGCAACCAGACCGTGGCCGATGGCAAGTTCTCCGAGGATGCTGATAGCGTGAGTGTGCATGGTGTGGAGCGCAAATTTCAACAATCGTCCTCGCAGAGTCTGCTGGCAACGTCACCGTGGAGCTCGAAGAGTTCCCGAACGTCCCGGGTTGCACCGACACCGGGTGTGTCACCAACGGGGTTGTCTAACCGGGGACGCAATCGATACAACACCTACAACGGGCCTGGTTCCGACTCGGCCACGTCGATTTCCACCGAGTCACACACCAGTCGACAGTCCCTCGTTACGCAGCATTCGGTGAGTGTAAGCTACGGTGTGTAGGCATTAGACCAGTGTAGCTTGAACGATACGTTGGTTCTGTCTTTTTAGCGATCTCGTACGCGAACAAGCACATTCCCGCCAAAGGTTTCATCGCAGGGCAGTCAGCAGGATAGGCAATCGGTGTCCCATTTTTTGTTGCAGCAACAGCAGCAGCAGCAGCCTCAGCAGCCGCAACAGCAATCATTCCAGCAGCAGCAACAAGCACCCCCAGATCTGCAACAGAGCCTGCAGCATTCGCGGCATCAACTGCAGCTACAGCAAAACCAGTTGCAGCAACAGTTACAGCAGCAACAACAACAACATCGTACGCAGCATCAGCAATTGCAACAATTGCAGGCTCACCAGCAACATCATCCTGCATCACAGCAACAACCTCAGCAGACTCAGCTCATGGTACAAACGATGTCGGAGCCGGATTTGTTTGGTCCTAAGATAGCG'},
 {'truncated': '', 'untruncated': ''},
"""

#______________________________________________________
# TEST 2: Check if '1312' will work as Release=1312 for genome access
#______________________________________________________
from cogent.db.ensembl import HostAccount, Genome, Species
account = HostAccount('localhost','vbuser','Savvas')
Species.amendSpecies('Anopheles albimanus','anopheles albimanus')
genome = Genome(Species='Anopheles albimanus',Release=73,account=account)



#______________________________________________________
# TEST 3: Check if UTR lengths according to locaiton obj agrees with UTR lengths according to actual str lenght of the utr sequence sampled
#______________________________________________________

import pickle
pickle_in = pickle.load(open('../data/sample_seqs/pickled/anopheles_gambiae_downstream.p','rb'))

utrLenghtOfObj_VS_utrLengthOfSeq = [(  (int(i.split('\t')[1].split(':')[2].split('-')[1]) - int(i.split('\t')[1].split(':')[2].split('-')[0])),int(i.split('\t')[2].replace('UtrLength:',''))) for i in pickle_in.keys()]

checkAgreement = [i for i in utrLenghtOfObj_VS_utrLengthOfSeq if not i[0]==i[1]]

if checkAgreement:
    print('Disagreements Exist: '+checkAgreement)
else:
    print('No Disagreements Exist')
    

#______________________________________________________
# TEST 4: old method only cared about up or downstream, new method is now sensitive to - and + strand, but is the new method in agreement with vectorbase genome browser?
#______________________________________________________

"""

...According to the...

    old method: 
        >AGAP004704 UtrCoord:2L:2253247-2253247:-1  UtrLength:0 UtrType:Synthetic
        ''

    new method:
        >AGAP004704 UtrCoord:2L:2253660-2253747:-1  UtrLength:87    UtrType:Synthetic
        'CTTTTTAAATACTATATTAGAAAAAGGTTAAGTTTATCTTATTATAATTTTTGTGATGCGTGTGGTACTCAATAAAACACGTGTGTT'

...presumably since...

    Overlapping Exons                  //---- -------- --------       
                                                           ||||
    NEW PERSPECTIVE                                        ------S------ + ---------G---------//   downstream +1  from the New Method's perspective
                                                |||||||||||||  <__87bp__>
    OLD PERSPECTIVE     //---------G--------- + ------S------                                      downstream +1  from the Old Method's perspective                        

                                               <0bp>

                            ...Thus 87bp remains under the new method, but 0bp remains under old method...

...So let's test if this 87bp is indeed a true sample sequence downstream of a -1 strand gene...

"""
# TESTING HERE USING VECTORBASE GENOME BROWSER & BLAST
"""

TEST VIA GENOME BROWSER:
     - I first enter the location of our NEW PERSPECTIVE sample location: '2L:2253660-2253747' into VectorBase's location search 
            - Location Search                         : https://www.vectorbase.org/Anopheles_gambiae/Location/View?db=core;g=AGAP004703;r=2L:2253660-2253747;t=AGAP004703-RA
     - Indeed we seems to be slap bang in between the gene in which we sampled from (AGAP004704) and some near-overlapping feature adjacent to our gene (AGAP004704)
     - By clicking on the gene in which we sampled from, it does indeed seem like we have sampled between the overlapping interface (link below) and the gene interface (link below)
            - Transcript of the overlapping interface : https://www.vectorbase.org/Anopheles_gambiae/Location/View?db=core;g=AGAP004703;r=2L:2248667-2253660;t=AGAP004703-RA      <= rightmost exon
            - Transcript of the gene interface        : https://www.vectorbase.org/Anopheles_gambiae/Transcript/Summary?db=core;g=AGAP004704;r=2L:2253748-2255495;t=AGAP004704-RA <= leftmost exon
            - Looking at both the overlap and gene    : https://www.vectorbase.org/Anopheles_gambiae/Gene/Summary?g=AGAP004704;r=2L:2253748-2255495;t=AGAP004704-RA

TEST VIA BLAST:
    - I BLAST the 87bp sampled sequence against An. gambiae:
        - BLAST RESULT                        : https://www.vectorbase.org/blast
            => There are no HSP matches that belong to a gene feature => good sign since it means we have truncated according to overlaps successfully        
            
    - So lets check the genome browser:
        - GENOME BROWSER                      : https://www.vectorbase.org/Anopheles_gambiae/Location/View?r=2L:2253631-2253777;time=1391619921
            => Good it seems we have landed slap bang in between the sampled-from-gene and an adjacent overlapping-gene

"""



