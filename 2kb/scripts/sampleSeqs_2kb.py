

# params:
upkb = 2000 # how much upstream seq to take 

import os

#______________________________________________#
# Register                           
#______________________________________________#
Release = 73
from cogent.db.ensembl import HostAccount
if 'ENSEMBL_ACCOUNT' in os.environ:
    host, username, password = os.environ['ENSEMBL_ACCOUNT'].split()
    account = HostAccount(host, username, password)
else:
    account = HostAccount('localhost','vbuser','Savvas')

#______________________________________________#
# Scan for genes to sample up/down stream seqs
#______________________________________________#
from cogent.db.ensembl import HostAccount, Genome, Species

Species.amendSpecies('Anopheles albimanus','anopheles albimanus')
Species.amendSpecies('Anopheles funestus','anopheles funestus')

ensemblGenome= Genome(Species='Anopheles albimanus', Release=Release, account=account)

# Define:
def scanGenes (ensemblGenome, geneId, featureType, scanFrom, scanTo):
    
    ## Example args:
    # ensemblGenome= Genome(Species=species, Release=Release, account=account)
    # scan_id     = 'ENSG00000012048'
    # scan_type   = 'gene'
    # scan_range  = (-100000,100000)

    print 'Scanning for genes...'
    scan_target = ensemblGenome.getGeneByStableId(StableId=geneId).Location
    scan_genes  = ensemblGenome.getFeatures(region=scan_target.resized(scanFrom,scanTo),feature_types=featureType)
    scan_geneIds= [i.StableId for i in scan_genes]
    return scan_geneIds
# Run:
geneIds = scanGenes(ensemblGenome,'AALB004251','gene',-100000000,100000000)
import pprint
#pprint.pprint([(i,j.Description) for i,j in enumerate(scan_genes)])

#______________________________________________#
# Sample up/down steam seqs from scanned genes
#______________________________________________#

import itertools

# Define:
def sampleSequences (sample_direction, sample_range, geneIds, sample_data={}):

    # ## Example Args:
    # sample_direction  = 'upstream' # 'downstream', 'both'
    # sample_range = 2000  # no. of bases to sample up/downstream
    # sample_data  = {} #Note: sampleSequences either starts from scratch (default) or builts upon previous data

    print 'Sampling sequences: '+str(sample_range)+'bp '+sample_direction+'...'

    from operator import attrgetter
    import pprint

    for geneId in geneIds:
        print '\t'+geneId
        gene = ensemblGenome.getGeneByStableId(StableId=geneId) # select gene
        geneLocation = gene.Location 
        sample_data[geneId] = {'geneLocation':geneLocation,'sampleLocation':{'untruncated':[],'truncated':[],'featuresIn':[]}}
     
        if sample_direction == 'upstream':
            sampleLocation = geneLocation.resized(-sample_range,-len(geneLocation)) #TODO: check if this is off by 1 or not
            overlap_features = list(ensemblGenome.getFeatures(region=sampleLocation,feature_types='gene'))
            #pprint.pprint(list(overlap_features))

            if not overlap_features: # skip to next loop if no features are found
                sample_data[geneId]['sampleLocation']['untruncated'] = sampleLocation 
                sample_data[geneId]['sampleLocation']['truncated']   = sampleLocation
                continue

            overlap_exons = [[[exon for exon in transcript.Exons] for transcript in transcripts] for transcripts in [feature.Transcripts for feature in overlap_features]] # exons Of Transcripts Of Features
            #pprint.pprint(overlap_exons)
            overlap_exons = list(itertools.chain(*overlap_exons)) # flattens
            #pprint.pprint(overlap_exons)
            overlap_exons = list(itertools.chain(*overlap_exons)) # flattens
            overlap_locations = sorted([exon.Location for exon in overlap_exons],key=attrgetter('End')) # 'End' <= 3' of feature vs. 5' of sample
            # TODO: ignore exons whose 5' end is > sample's 3' end (not technically overlapping)

            # DEBUG BREAKPOINT # If something went wrong.. this can be used to escape the function and return desired vars for debug
            try:
                overlap_limit     = overlap_locations[-1] # the last overlapping feature beyond which there are no overlaps w/i sample
            except IndexError:
                print('Mission Abort!!!...Activating Debug route...')
                return overlap_exons, overlap_features

            if overlap_limit.End > sampleLocation.End: # sometimes the entire sample is nested within an overlapping feature, then we have no sample
                print('\t\tNo free upstream sequence available for: '+geneId)
                sampleLocation_truncated= None
            else:
                sampleLocation_truncated= sampleLocation.resized(+(overlap_limit.End-sampleLocation.Start),0) # Truncate the sampled upstream sequence (promoter to prevent overlap)

            sample_data[geneId]['sampleLocation']['untruncated'] = sampleLocation 
            sample_data[geneId]['sampleLocation']['truncated']   = sampleLocation_truncated 
            sample_data[geneId]['sampleLocation']['featuresIn']  = overlap_features

        elif sample_direction == 'downstream':
            print '\t'+sample_direction+'...'
            sampleLocation      = geneLocation.resized(len(geneLocation),+sample_range)
            overlap_features    = [feature for feature in ensemblGenome.getFeatures(region=sampleLocation,feature_types='gene')]
            overlap_locations   = sorted([feature.Location for feature in overlap_features],key=attrgetter('Start'))    # 'Start' <= 5' of feature vs. 3' of sample
    return sample_data, feature
# Run:

sample_data, overlap_exons = sampleSequences('upstream', 2000, geneIds)
#overlap_exons, overlap_features = sampleSequences('upstream', 2000, geneIds)

import pprint
pprint.pprint(sample_data)