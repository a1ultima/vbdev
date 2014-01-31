# IMPORT
from operator import attrgetter
import itertools


# DEFINE FUNCTIONS

# #______________________________________________#
# # Scan for genes to sample up/down stream seqs
# #______________________________________________#
# def scanGenes (ensembl_genome, geneId, featureType, scanFrom, scanTo):
#     """ Example args:
#     ensembl_genome= Genome(Species=species, Release=Release, account=account)
#     scan_id     = 'ENSG00000012048'
#     scan_type   = 'gene'
#     scan_range  = (-100000,100000) 
#     """
#     print 'Scanning for genes...'
#     scan_target = ensembl_genome.getGeneByStableId(StableId=geneId).Location
#     scan_genes  = ensembl_genome.getFeatures(region=scan_target.resized(scanFrom,scanTo),feature_types=featureType)
#     scan_geneIds= [i.StableId for i in scan_genes]
#     return scan_geneIds


#______________________________________________#
# Sample up/down steam seqs from scanned genes
#______________________________________________#
# def scanGenes (ensembl_genome, geneId, featureType, scanFrom, scanTo):
#     """ Example args:
#     ensembl_genome= Genome(Species=species, Release=Release, account=account)
#     scan_id     = 'ENSG00000012048'
#     scan_type   = 'gene'
#     scan_range  = (-100000,100000) 
#     """
#     print 'Scanning for genes...'
#     scan_target = ensembl_genome.getGeneByStableId(StableId=geneId).Location
#     scan_genes  = ensembl_genome.getFeatures(region=scan_target.resized(scanFrom,scanTo),feature_types=featureType)
#     scan_geneIds= [i.StableId for i in scan_genes]
#     return scan_geneIds

#______________________________________________#
# Genome Setup
#______________________________________________#

def setupGenome( species, db_host=None, db_user=None, db_pass=None, db_release=None ):
    """ Example Args:
    species     = 'Anopheles gambiae'
    db_host     = 'localhost'
    db_user     = 'vbuser'
    db_pass     = 'Savvas'
    db_release  = 73
    """
    if not db_host:
        account=None
    if not db_release:
        db_release=73
    from cogent.db.ensembl import HostAccount, Genome, Species
    account = HostAccount(db_host,db_user,db_pass)
    Species.amendSpecies(species,species)
    genome = Genome(Species=species,Release=db_release,account=account)
    return genome

#______________________________________________#
# Sample up/down steam seqs from scanned genes
#______________________________________________#
def sampleSequences_read ( ensembl_genome, sample_direction='upstream', sample_range=2000, sample_data={}, sample_seqs={}):
    """ Example Args:
    sample_direction  = 'upstream' # 'downstream', 'both'
    sample_range = 2000  # no. of bases to sample up/downstream
    sample_data  = {} #Note: sampleSequences either starts from scratch (default) or builts upon previous data 
    """
    print 'Sampling sequences: '+str(sample_range)+'bp '+sample_direction+'...'
    genes   = ensembl_genome.getGenesMatching(BioType='protein_coding')   # queries ensembl genome for all protein coding genes
    geneIds = [gene.StableId for gene in genes]  # grab all gene ids
    geneIds = geneIds[0:100]    # ANDY: testing for just 100

    for geneId in geneIds:
        print '\t'+geneId
        geneLocation                = ensembl_genome.getGeneByStableId(StableId=geneId).Location # select gene
        #geneLocation        = gene.Location 
        sample_data[geneId] = {'geneLocation':geneLocation,'sampleLocation':{'untruncated':[],'truncated':[],'featuresIn':[]}} # initiate storage
        sample_seqs[geneId] = {'untruncated':str(),'truncated':str()}

        if sample_direction=='upstream':
            sampleLocation  = geneLocation.resized(-sample_range,-len(geneLocation)) # focus on the sampled region, #TODO: check if off by 1 @@@@@@@@@@@
        elif sample_direction=='downstream':
            sampleLocation  = geneLocation.resized(+len(geneLocation),+sample_range) 
        overlap_features    = list(ensembl_genome.getFeatures(region=sampleLocation,feature_types='gene'))

        if not overlap_features:            # keep full sample & skip to next loop if no features are found
            # a) NO OVERLAPS PROCEDURE:
            print('\t\tNo overlapping features')
            sample_data[geneId]['sampleLocation']['untruncated']= sampleLocation  # Alterrnate route 1 --v
            sample_data[geneId]['sampleLocation']['truncated']  = sampleLocation
            sample_seqs[geneId]['untruncated']                  = str(ensembl_genome.getRegion(sampleLocation).Seq)
            sample_seqs[geneId]['truncated']                    = str(ensembl_genome.getRegion(sampleLocation).Seq)
        else:
            # b) OVERLAPS PROCEDURE:
            overlap_exons       = [[[exon for exon in transcript.Exons] for transcript in transcripts] for transcripts in [feature.Transcripts for feature in overlap_features]] # exons Of Transcripts Of Features
            overlap_exons       = list(itertools.chain(*list(itertools.chain(*overlap_exons)))) # flattens from: 3d->2d->1d
            overlap_locations   = [exon.Location for exon in overlap_exons]

            # TRUNCATION STEP 1:    remove exons that are so far down/upstream of sample they are beyond sample overlap // WARN: can be redundant?
            #
            #   overlapping exons:      5'----3'    5'----3'         5'----3'       5'----3' 
            #                              ||||      ||||||            xxxx           xxxx  <--notice we dont want the xxx's to truncate away our Sample
            #   sample:                    5'-------S---------3'+5'-----------G------------------//  S = sample, G = gene of sample
            #
            if sample_direction=='upstream':
                overlap_locations = [i for i in overlap_locations if i.Start < sampleLocation.End]
            elif sample_direction=='downstream':
                overlap_locations = [i for i in overlap_locations if i.End > sampleLocation.Start]

            # TRUNCATION STEP 2:    care only about the furthest downstream exon who overlaps w/ sample
            #
            #   overlapping exons:     >> 5'---3' <<   5'---3' 5'---3' 5'---3' 5'---3'
            #
            if sample_direction=='upstream':
                overlap_locations = sorted(overlap_locations,key=attrgetter('End')) # Sort by 'End' <= 3' of feature vs. 5' of sample
                overlap_limit = overlap_locations[-1] # limiting overlapping feature beyond which there are no overlaps w/i sample
            elif sample_direction=='downstream':
                overlap_locations = sorted(overlap_locations,key=attrgetter('Start')) # 'End' <= 3' of feature vs. 5' of sample
                overlap_limit = overlap_locations[0]  

            # TRUNCATION STEP 3:    
            #
            #   overlapping exon: 5'---------------------3'
            #                         |||||||||||||||||||  <--notice the entire sample must be removed 
            #   sample:               5'----S----3'+5'--------G--------//  
            #
            if sample_direction == 'upstream':
                if overlap_limit.End > sampleLocation.End:          # if overlapping feature spans the entire sample
                    print('\t\tNo free upstream sequence available for: '+geneId)
                    sampleLocation_truncated= None
                else:
                    sampleLocation_truncated= sampleLocation.resized(+(overlap_limit.End-sampleLocation.Start),0) # Truncate the sampled upstream sequence (promoter to prevent overlap)
            elif sample_direction == 'downstream':
                if overlap_limit.Start < sampleLocation.Start:      # if overlapping feature spans the entire sample
                    print('\t\tNo free upstream sequence available for: '+geneId)
                    sampleLocation_truncated= None
                else:
                    sampleLocation_truncated= sampleLocation.resized(0,-(sampleLocation.End - overlap_limit.Start)) # Truncate the sampled upstream sequence (promoter to prevent overlap)

            # STORE DATA 
            sample_data[geneId]['sampleLocation']['untruncated']= sampleLocation # Alternate route 2   --^
            sample_data[geneId]['sampleLocation']['truncated']  = sampleLocation_truncated 
            sample_data[geneId]['sampleLocation']['featuresIn'] = overlap_features
            sample_seqs[geneId]['untruncated']                  = str(ensembl_genome.getRegion(sampleLocation).Seq)
            if sampleLocation_truncated:
                sample_seqs[geneId]['truncated']                = str(ensembl_genome.getRegion(sampleLocation_truncated).Seq)
            else:
                sample_seqs[geneId]['truncated']                = None

    #return sample_data, sample_seqs    #ANDY_01/27
    return sample_seqs  # yields the dictionary

#______________________________________________#
# write sequences into either a pickle or fasta
#______________________________________________#
def sampleSequences_write(ensembl_genome, sample_data, fasta_it=True, pickle_it=True):

    samples  = {}
    geneIds  = sample_data.keys()
    species  = ensembl_genome.Species.replace(' ','_').lower()

    for geneId in geneIds:
        seq         = sample_data[geneId]['truncated']
        if seq:
            seq_length = len(seq)
        else:
            seq_length = 0
            seq = ''
        header = geneId+'_'+str(seq_length)+'bp'
        samples[header] = seq

    if pickle_it:
        import pickle
        pickle_path = '../data/sample_seqs/pickled/'
        pickle.dump(samples,open(pickle_path+species+'.p','wb'))

    if fasta_it:
        from cogent import LoadSeqs, DNA
        fasta_path      = '../data/sample_seqs/fasta/'
        fasta_file      = open(fasta_path+species+'.fasta','wb')

        #return samples
        samples_fasta   = LoadSeqs(data=samples,moltype=DNA,aligned=False).toFasta()

        fasta_file.write(samples_fasta)
        fasta_file.close()

samples = sampleSequences_write(genome,sample_data)

#-------------------------------------------------------------------------------------------------------------------------------------
# OLD METHODS

#______________________________________________#
# Buffer the sampling from a list of genes
#______________________________________________#
# def sampleSeqs_batch (...):

# EXECUTE
# need an iteratively naming system that works as the gene scanner iterates upwards 
# account         = ('localhost','vbuser','Savvas')
# species         = 'Anopheles albimanus'
# release         = 73

# #starterGene     = 'AALB004251'
# sample_range    = 2000
# sample_direction= 'upstream'
# account         = HostAccount(account[0],account[1],account[2]) 
# Species.amendSpecies(species,species)
# ensembl_genome   = Genome(Species=species, Release=release, account=account)
# geneList        = scanGenes(ensembl_genome,starterGene,'gene',-99999999999,99999999999) # Grabs a list of genes to sample from --> ideally if we just had all the genes in a list...
# print geneList


# BATCH IT METHOD... but do we need?
#   This buffers the procedure to X genes saved to pickle each loop
# print('Buffer Mode: ON')
# bufferSize          = 100 # No. of genes at a time     to sample from
# geneList_chunked    = list(chunks(geneList,bufferSize))
# pathout             = '../data/'
# nGenes_from         = 1
# for i,chunk in enumerate(geneList_chunked):
#     nGenes_to   = nGenes_from + len(chunk)
#     fileout     = 'batch'+str(i)+'_genes_'+str(nGenes_from)+'-to-'+str(nGenes_to)
#     nGenes_from = nGenes_to
#     print(fileout)
#     sample_data, sample_seqs = sampleSequences(sample_direction, sample_range, chunk)   #ANDY_01/27
#     pickle.dump( sample_seqs, open( pathout+fileout+".p", "wb" ) ) # Save a dictionary into a pickle file.


