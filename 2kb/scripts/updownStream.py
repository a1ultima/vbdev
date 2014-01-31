# IMPORT
from operator import attrgetter
import itertools

import difflib


# DEFINE FUNCTIONS

#______________________________________________#
# Genome Setup
#______________________________________________#

def setupGenome( species, db_host=None, db_user=None, db_pass=None, db_release=None ):
    """ 
    NOTES:
        setup the ensembl_genome object using pycogent, this object has methods described here: http://pycogent.org/examples/query_ensembl.html
    
    ARGUMENTS:
        species     = 'Anopheles gambiae'   #string: this needs to match the mysql databases in vbdev: >mysql -hlocalhost -uvbuser -pSavvas
        db_host     = 'localhost'           #string: ''
        db_user     = 'vbuser'              #string: ''
        db_pass     = 'Savvas'              #string: ''
        db_release  = 73                    #integer: the realease versions we use are all 73
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
def sampleSequences_read( ensembl_genome, sample_direction='upstream', sample_range=2000, sample_data={}, sample_seqs={}):
    """ 
    Notes:
        takes an ensemble_genome via setupGenome() and starts sampling up/downstream sequences from all the genes within

    Args:
        ensembl_genome = setupGenome( ... )                             # ensembl_genome object:  generated by setupGenome( ... ) 
        sample_direction  = 'upstream' # 'downstream', 'both'           # string:       sample Xbp upstream or downstream from each gene??
        sample_range = 2000  # no. of bases to sample up/downstream     # int:          how many bp of sequence to sample?
        sample_data  = {}                                               # dictionary:   starts from scratch (default) or recurse upon previous data 
    """
    print 'Sampling sequences: '+str(sample_range)+'bp '+sample_direction+'...'
    genes   = ensembl_genome.getGenesMatching(BioType='protein_coding')   # queries ensembl genome for all protein coding genes
    geneIds = [gene.StableId for gene in genes]  # grab all gene ids
    geneIds = geneIds[0:10]    # ANDY: testing for just 100

    for count,geneId in enumerate(geneIds):
        print '\t'+str(count)+' '+geneId
        gene        = ensembl_genome.getGeneByStableId(StableId=geneId) # select gene
        geneLocation= gene.Location # coordinates for whole gene

        # Assuming all genes have utr3 annotated:       <- WARN: not always true?

        #gene_exons  = [[e.Location for e in t.Exons] for t in gene.Transcripts] # coordinates for limiting exon fo the gene

    # EXON TARGET               # Sampling starts from the limiting exon of the current gene -> and spans according to sample_range
        if sample_direction=='upstream':
            
            print('cool')

            # gene_exons          = sorted(list(itertools.chain(*gene_exons)),key=attrgetter('Start')) # flattens from: 2d->1d
            # gene_limiting_exon  = gene_exons[0]     # exon closest to the 5' end of sample
            # sampleLocation      = gene_limiting_exon.resized(-sample_range-1,-len(gene_limiting_exon)-1) # focus on the sampled region, #TODO: check if off by 1 @@@@@@@@@@@
        elif sample_direction=='downstream':

            seq_transcript = str(gene.getLongestCdsTranscript().Seq)
            seq_utr        = str(gene.getLongestCdsTranscript().Utr3)

            matcher = difflib.SequenceMatcher(a=seq_transcript, b=seq_utr) # Matches coordinates of a string and its substring
            match   = matcher.find_longest_match(0, len(matcher.a), 0, len(matcher.b))

            location_sample     = gene.getLongestCdsTranscript().Location.resized(  match.a,   # <- where the utr begins amtching                         
                                                                                    match.b+(sample_range-match.size)) # <- where the utr needs to extend for specified sample_range
            location_transcript = gene.getLongestCdsTranscript().Location.resized(0,match.a)   # <- truncate the utr away

            seq_transcript  = str(ensembl_genome.getRegion(location_transcript).Seq)
            seq_sample      = str(ensembl_genome.getRegion(location_sample).Seq)

            sampleLocation = location_sample
            # gene_exons          = sorted(list(itertools.chain(*gene_exons)),key=attrgetter('End'))  # flattens from: 2d->1d
            # gene_limiting_exon  = gene_exons[-1]    # exon closest to 3' end of sample
            # sampleLocation      = gene_limiting_exon.resized(+len(gene_limiting_exon)+1,+sample_range+1) 

    # INITIATE STORAGE:         ANDY: concatenate just the limiting exon sequence to the sampleseq?
        #sample_seqs[geneId] = {'untruncated':str(),'truncated':str(),'gene_seq':str(ensembl_genome.getRegion(gene_limiting_exon).Seq),'sample_location':str(),'gene_location':':'.join(str(gene_limiting_exon).split(':')[2:])}
        sample_seqs[geneId] = {'untruncated':str(),'truncated':str(),'gene_seq':str(seq_transcript),'sample_location':str(),'gene_location':':'.join(str(location_transcript).split(':')[2:])}

    # FEATURES OVERLAP WITH SAMPLE? 
        overlaps_flag           = None
        overlap_features        = list(ensembl_genome.getFeatures(region=sampleLocation,feature_types='gene'))
        if overlap_features:
            overlaps_flag       = True
            overlap_exons       = [[[exon for exon in transcript.Exons] for transcript in transcripts] for transcripts in [feature.Transcripts for feature in overlap_features]] # exons Of Transcripts Of Features
            overlap_exons       = list(itertools.chain(*list(itertools.chain(*overlap_exons)))) # flattens from: 3d->2d->1d
            overlap_locations   = [exon.Location for exon in overlap_exons]
            # PRUNE FALSE OVERLAPS:    remove exons that are so far down/upstream of sample they are beyond sample overlap // WARN: can be redundant?
            #
            #   overlapping exons:      5'----3'    5'----3'         5'----3'       5'----3' 
            #                              ||||      ||||||            xxxx           xxxx  <- notice we dont want the xxx's to truncate away our Sample
            #   sample:                    5'-------S---------3'+5'-----------G------------------//  S = sample, G = gene of sample
            overlap_locations   = [i for i in overlap_locations if (i.Start < sampleLocation.End) and (i.End > sampleLocation.Start)] # remove "too-far-downstream" features
            if overlap_locations==[]:
                overlaps_flag = False
        else:
            overlaps_flag = False

    # SAMPLE TRUNCATION PROCEDURES
        if not overlaps_flag:            # keep full sample & skip to next loop if no features are found
            # a) NO OVERLAPS? -> Store whole sample

            # print('\t\tOvelapping Features: NO')
            sample_seqs[geneId]['untruncated']      = sample_seqs[geneId]['truncated'] = str(ensembl_genome.getRegion(sampleLocation).Seq)
            sample_seqs[geneId]['sample_location']  = ':'.join(str(sampleLocation).split(':')[2:]) # genome coordinates summary
        else:
            # b) OVERLAPS? -> Truncate sample accordingly
            
            # print('\t\tOverLapping Features: YES')

            # TRUNCATION STEP 2:    care only about the furthest downstream exon who overlaps w/ sample
            #
            #   overlapping exons:  Only want this one!->  >>>5'---3'<<<   5'---3' 5'---3' 5'---3' 5'---3'
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
                    #print('\t\tComplete overlap!')
                    sampleLocation_truncated= None
                else:
                    sampleLocation_truncated= sampleLocation.resized(+(overlap_limit.End-sampleLocation.Start),0) # Truncate the sampled upstream sequence (promoter to prevent overlap)
            elif sample_direction == 'downstream':
                if overlap_limit.Start < sampleLocation.Start:      # if overlapping feature spans the entire sample
                    #print('\t\tComplete overlap!')
                    sampleLocation_truncated= None
                else:
                    sampleLocation_truncated= sampleLocation.resized(0,-(sampleLocation.End - overlap_limit.Start)) # Truncate the sampled upstream sequence (promoter to prevent overlap)
        # STORE DATA...
            sample_seqs[geneId]['sample_location']              = ':'.join(str(sampleLocation).split(':')[2:]) # genome coordinates summary
            sample_seqs[geneId]['untruncated']                  = str(ensembl_genome.getRegion(sampleLocation).Seq)
            if sampleLocation_truncated:
                sample_seqs[geneId]['truncated']                = str(ensembl_genome.getRegion(sampleLocation_truncated).Seq)
            else:
                sample_seqs[geneId]['truncated']                = ''

    #return sample_data, sample_seqs    #ANDY_01/27
    return sample_seqs, sample_direction  # yields the dictionary

#______________________________________________#
# write sequences into either a pickle or fasta
#______________________________________________#
def sampleSequences_write(ensembl_genome, sample_read, fasta_it=True, pickle_it=True, include_genes=False):
    """ 
    Notes:
        Stores into files (either .fasta or .p) the sampled sequence data that was generated into a dict by sampleSequences_read()

    Args:
        ensembl_genome  = setupGenome( ... )                # ensembl genome object: http://pycogent.org/examples/query_ensembl.html
        sample_read     = sampleSequences_read( ... )       # tuple: ( dict of sampled seqs, sample_direction e.g. 'upstream')
        fasta_it        = True                              # bool: generate a fasta file for all seqs? 
        pickle_it       = False                             # bool: store a python dict with all the sequences as file?
        include_genes   = False                             # bool: concatenate sampled seqs onto gene seqs?
    """
    samples         = {}
    sample_data     = sample_read[0]
    sample_direction= sample_read[1]
    geneIds         = sample_data.keys()
    species         = ensembl_genome.Species.replace(' ','_').lower()
    for geneId in geneIds:
        if include_genes:   # concatenate sample_seq to gene_seq?
            gene_seq        = sample_data[geneId]['gene_seq']
            sample_seq      = sample_data[geneId]['truncated']
            if sample_direction == 'downstream':
                total_seq       = gene_seq+sample_seq 
            elif sample_direction == 'upstream':
                total_seq       = sample_seq+gene_seq
            gene_length     = len(gene_seq)
            sample_length   = len(sample_seq)
            header = geneId+'\tGeneCoord:'+sample_data[geneId]['gene_location']+'\tUTRCoord:'+sample_data[geneId]['sample_location']+'\tGeneLength:'+str(gene_length)+'\tUTRLength:'+str(sample_length)
        else:
            total_seq       = sample_data[geneId]['truncated']
            total_length    = len(total_seq)
            header = geneId+'\tUTRCoord:'+sample_data[geneId]['sample_location']+'\tUTRLength:'+str(total_length)
        samples[header] = total_seq
    if pickle_it:
        import pickle
        pickle_path     = '../data/sample_seqs/pickled/'
        pickle.dump(samples,open(pickle_path+species+'_'+sample_direction+'.p','wb'))
    if fasta_it:
        from cogent import LoadSeqs, DNA
        fasta_path      = '../data/sample_seqs/fasta/'
        fasta_file      = open(fasta_path+species+'_'+sample_direction+'.fasta','wb')
        samples_fasta   = LoadSeqs(data=samples,moltype=DNA,aligned=False).toFasta()
        fasta_file.write(samples_fasta)
        fasta_file.close()
# SUMMARY STATS
    if include_genes:
        headers = [i.split('\t') for i in samples.keys()]
        samplelength_distribution = [int(i[4].replace('UTRLength:','')) for i in headers]
    else:
        headers = [i.split('\t') for i in samples.keys()]
        samplelength_distribution = [int(i[2].replace('UTRLength:','')) for i in headers]
    print('Complete!')
    return headers, samplelength_distribution

#______________________________________________#
# samples from a list of genomes
#______________________________________________#
def sampleAllSpecies(species_list,sample_directions):
    """
    Notes:
        samples sequences up/downstream from each gene of each genome of a list of species generating fasta files
    Args:
        species_list        = ['Anopheles gambiae', 'Aedes aegypti']    #list: species names to sample from
        sample_directions   = ['upstream','downstream']                 #list: list of directions to sample seqs from
    """
    for sample_direction in sample_directions:
        for species in species_list:
            print(species)
            genome          = setupGenome(              species, db_host='localhost', db_user='vbuser', db_pass='Savvas', db_release=73 )
            samples_read    = sampleSequences_read (    genome, sample_direction=sample_direction, sample_range=200, sample_data={}, sample_seqs={})
            samples_write   = sampleSequences_write(    genome, samples_read, fasta_it=True, pickle_it=True, include_genes=True)
            import pprint
            pprint.pprint(samples_write[1]) # show the header

#-------------------------------------------------------------------------------------------------------------------------------------
# RUN
#-------------------------------------------------------------------------------------------------------------------------------------

species_list = [    'Aedes aegypti',
                'Anopheles albimanus',
                'Anopheles arabiensis',
                'Anopheles atroparvus',
                'Anopheles christyi',
                'Anopheles culicifacies',
                'Anopheles darlingi',
                'Anopheles dirus',
                'Anopheles epiroticus',
                'Anopheles farauti',
                'Anopheles funestus',
                'Anopheles gambiae',
                'Anopheles maculatus',
                'Anopheles melas',
                'Anopheles merus',
                'Anopheles minimus',
                'Anopheles quadriannulatus',
                'Anopheles sinensis'
                #,'Anopheles stephensiI',           # this is broken..?
                #'Anopheles stephensi'              # ''
                ]
sample_directions = ['downstream','upstream']

#sampleAllSpecies(species_list,sample_directions)

#-------------------------------------------------------------------------------------------------------------------------------------
# TESTING
#-------------------------------------------------------------------------------------------------------------------------------------

genome = setupGenome(              'Anopheles gambiae', db_host='localhost', db_user='vbuser', db_pass='Savvas', db_release=73 )
gene   = sampleSequences_read (    genome, sample_direction='downstream', sample_range=200, sample_data={}, sample_seqs={})


#-------------------------------------------------------------------------------------------------------------------------------------
# OLD METHODS
#-------------------------------------------------------------------------------------------------------------------------------------

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
