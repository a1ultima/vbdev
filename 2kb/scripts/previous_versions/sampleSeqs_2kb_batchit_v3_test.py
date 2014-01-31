# IMPORT
import os
import pickle
from cogent.db.ensembl import HostAccount, Genome, Species
from cogent.db.ensembl.util import NoItemError
from operator import attrgetter
import pprint
from cogent.db.ensembl import HostAccount
import itertools

# DEFINE FUNCTIONS

#______________________________________________#
# Scan for genes to sample up/down stream seqs
#______________________________________________#
def scanGenes (ensemblGenome, geneId, featureType, scanFrom, scanTo):
    """ Example args:
    ensemblGenome= Genome(Species=species, Release=Release, account=account)
    scan_id     = 'ENSG00000012048'
    scan_type   = 'gene'
    scan_range  = (-100000,100000) 
    """
    
    print 'Scanning for genes...'

    scan_target = ensemblGenome.getGeneByStableId(StableId=geneId).Location
    scan_genes  = ensemblGenome.getFeatures(region=scan_target.resized(scanFrom,scanTo),feature_types=featureType)
    scan_geneIds= [i.StableId for i in scan_genes]
    
    return scan_geneIds

#______________________________________________#
# Sample up/down steam seqs from scanned genes
#______________________________________________#
def sampleSequences (sample_direction, sample_range, geneIds, sample_data={}, sample_seqs={}):
    """ Example Args:
    sample_direction  = 'upstream' # 'downstream', 'both'
    sample_range = 2000  # no. of bases to sample up/downstream
    sample_data  = {} #Note: sampleSequences either starts from scratch (default) or builts upon previous data 
    """

    print 'Sampling sequences: '+str(sample_range)+'bp '+sample_direction+'...'

    print(geneIds)

    for geneId in geneIds:
        print '\t'+geneId
        gene = ensemblGenome.getGeneByStableId(StableId=geneId) # select gene
        geneLocation = gene.Location 
        sample_data[geneId] = {'geneLocation':geneLocation,'sampleLocation':{'untruncated':[],'truncated':[],'featuresIn':[]}}
        sample_seqs[geneId] = {'untruncated':str(),'truncated':str()}
     
        if sample_direction == 'upstream':
            sampleLocation = geneLocation.resized(-sample_range,-len(geneLocation)) #TODO: check if this is off by 1 or not
            overlap_features = list(ensemblGenome.getFeatures(region=sampleLocation,feature_types='gene'))

            if not overlap_features:            # keep full sample & skip to next loop if no features are found
                # a) NO OVERLAPS PROCEDURE:     # TODO: print msg that there are no overlaps
                sample_data[geneId]['sampleLocation']['untruncated']= sampleLocation  # Alterrnate route 1 --v
                sample_data[geneId]['sampleLocation']['truncated']  = sampleLocation
                sample_seqs[geneId]['untruncated']                  = str(ensemblGenome.getRegion(sampleLocation).Seq)
                sample_seqs[geneId]['truncated']                    = str(ensemblGenome.getRegion(sampleLocation).Seq)
            else:
                # b) OVERLAPS PROCEDURE:
                overlap_exons = [[[exon for exon in transcript.Exons] for transcript in transcripts] for transcripts in [feature.Transcripts for feature in overlap_features]] # exons Of Transcripts Of Features
                overlap_exons = list(itertools.chain(*overlap_exons)) # flattens
                overlap_exons = list(itertools.chain(*overlap_exons)) # flattens
                overlap_locations = sorted([exon.Location for exon in overlap_exons],key=attrgetter('End')) # 'End' <= 3' of feature vs. 5' of sample
                # TODO: ignore exons whose 5' end is > sample's 3' end (not technically overlapping)
                
                try:    # DEBUG BREAKPOINT                      # If something went wrong.. this can be used to escape the function and return desired vars for debug
                    overlap_limit     = overlap_locations[-1]   # the last overlapping feature beyond which there are no overlaps w/i sample
                except IndexError:
                    print('Mission Abort!!!...Activating Debug route...')
                    return overlap_exons, overlap_features

                #   overlapping exon 5'---------------------3'
                #                           ||||||||||||||||   <-- notice the entire sample must be removed 
                #   sample                  5'----3'-5'---------------3'
                #                              ^upstream sample   ^gene of sample
                if overlap_limit.End > sampleLocation.End:      # sometimes the entire sample is nested within an overlapping feature, then we have no sample
                    print('\t\tNo free upstream sequence available for: '+geneId)
                    sampleLocation_truncated= None
                else:
                    sampleLocation_truncated= sampleLocation.resized(+(overlap_limit.End-sampleLocation.Start),0) # Truncate the sampled upstream sequence (promoter to prevent overlap)
                # Stores the data away...
                sample_data[geneId]['sampleLocation']['untruncated']= sampleLocation # Alternate route 2   --^
                sample_data[geneId]['sampleLocation']['truncated']  = sampleLocation_truncated 
                sample_data[geneId]['sampleLocation']['featuresIn'] = overlap_features
                sample_seqs[geneId]['untruncated']                  = str(ensemblGenome.getRegion(sampleLocation).Seq)
                if sampleLocation_truncated:
                    sample_seqs[geneId]['truncated']                = str(ensemblGenome.getRegion(sampleLocation_truncated).Seq)
                else:
                    sample_seqs[geneId]['truncated']                = None
        elif sample_direction == 'downstream':
            print '\t'+sample_direction+'...'
            sampleLocation      = geneLocation.resized(len(geneLocation),+sample_range)
            overlap_features    = [feature for feature in ensemblGenome.getFeatures(region=sampleLocation,feature_types='gene')]
            overlap_locations   = sorted([feature.Location for feature in overlap_features],key=attrgetter('Start'))    # 'Start' <= 5' of feature vs. 3' of sample

    return sample_data, sample_seqs    #ANDY_01/27
    #return sample_data, sample_seqs

#______________________________________________#
# Buffer the sampling from a list of genes
#______________________________________________#
# def sampleSeqs_batch (...):

# EXECUTE
# need an iteratively naming system that works as the gene scanner iterates upwards 
account         = ('localhost','vbuser','Savvas')
species         = 'Anopheles albimanus'
release         = 73
starterGene     = 'AALB004251'
sample_range    = 2000
sample_direction= 'upstream'

account         = HostAccount(account[0],account[1],account[2]) 
Species.amendSpecies(species,species)
ensemblGenome   = Genome(Species=species, Release=release, account=account)
geneList        = scanGenes(ensemblGenome,starterGene,'gene',-99999999999,99999999999) # Grabs a list of genes to sample from --> ideally if we just had all the genes in a list...

print geneList

#pprint.pprint([(i,j.Description) for i,j in enumerate(scan_genes)])

def chunks(l, n):
    """ Yield successive n-sized chunks from l. """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

print('Buffer Mode: ON')
bufferSize          = 100 # No. of genes at a time     to sample from
geneList_chunked    = list(chunks(geneList,bufferSize))
pathout             = '../data/'
nGenes_from         = 1

for i,chunk in enumerate(geneList_chunked):
    nGenes_to   = nGenes_from + len(chunk)
    fileout     = 'batch'+str(i)+'_genes_'+str(nGenes_from)+'-to-'+str(nGenes_to)
    nGenes_from = nGenes_to
    print(fileout)
    #bug = sampleSequences(sample_direction, sample_range, chunk)   #ANDY_01/27
    #break
    sample_data, sample_seqs = sampleSequences(sample_direction, sample_range, chunk)   #ANDY_01/27
    pickle.dump( sample_seqs, open( pathout+fileout+".p", "wb" ) ) # Save a dictionary into a pickle file.


