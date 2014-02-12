# IMPORT
from operator import attrgetter
import itertools

#import difflib
import re


# TODOs -->                                                                                                                                 # v-- TODOs align HERE --v

# DEFINE FUNCTIONS

#______________________________________________#
# Genome Setup
#______________________________________________#
def setupGenome( species, db_host=None, db_user=None, db_pass=None, db_release=None ):
    """ 
    Notes:
        setup the ensembl_genome object using pycogent, this object has methods described here: http://pycogent.org/examples/query_ensembl.html
    
    Args:
        species     = 'Anopheles gambiae'   #string: this needs to match the mysql databases in vbdev: >mysql -hlocalhost -uvbuser -pSavvas
        db_host     = 'localhost'           #string: ^^
        db_user     = 'vbuser'              #string: ^^
        db_pass     = 'Savvas'              #string: ^^
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
def sampleSequences_read( ensembl_genome, sample_direction='upstream', take_annotated_utr=False, sample_range=500, masking = 'none', sample_data={}, sample_seqs={}, test_it=False):
    """ 
    Notes:
        takes an ensemble_genome via setupGenome() and starts sampling up/downstream sequences from all the genes within

    Args:
        ensembl_genome      = setupGenome( ... )                            # ensembl_genome object:  generated by setupGenome( ... ) 
        sample_direction    = 'upstream' # 'downstream', 'both'         	# string:       sample Xbp upstream or downstream from each gene??
        take_annotated_utr  = False                                         # bool:         do we accept already-annotated utr sequences of len>3bp, rather than syntehsizing it?
        sample_range        = 2000  # no. of bases to sample up/downstream  # int:          how many bp of sequence to sample?
        masking             = 'soft'                                        # string:       'soft' / 'hard' / 'none'  what kind of masking procedure do we use?
        sample_data         = {}                                            # dictionary:   starts from scratch (default) or recurse upon previous data 
    """

    print 'Sampling sequences: '+str(sample_range)+'bp '+sample_direction+'...'
    sample_failures = []
    genes           = ensembl_genome.getGenesMatching(BioType='protein_coding')   # queries ensembl genome for all protein coding genes
    geneIds         = [gene.StableId for gene in genes]  # grab all gene ids
    
    annotated_sample_flag = None
    
    if test_it == True:
        geneIds = geneIds[0:10]    # ANDY: testing for just 100
        #geneIds = ['AGAP000002']

    for count,geneId in enumerate(geneIds):
        if not test_it:
            if count%10==0:
                print '\t'+str(count)+' '+geneId
        else:
            print '\t'+str(count)+' '+geneId
        gene = ensembl_genome.getGeneByStableId(StableId=geneId) # select gene
        #geneLocation= gene.Location    # coordinates for whole gene

    # EXON TARGET                       # Sampling starts from the limiting exon of the current gene  and spans according to sample_range
        if sample_direction=='upstream':

            transcript          = gene.getLongestCdsTranscript()
            seq_transcript      = str(transcript.Seq)
            location_transcript = transcript.Location
            sample_type         = 'Synthetic'                   # for upstream sampling, we purely synthesize the data - no care for annotated

            if test_it:
                print('\t\tUTR Synthetic')

        # ORIENTATION?
            # (-) strand
            if location_transcript.Strand == -1:
                if test_it:
                    print('\t\t-')
                location_sample = location_transcript.resized( 0+len(seq_transcript), 0+sample_range )  # (shift to transcript_end, shift to sampling range)
            # (+) strand
            else:
                if test_it:
                    print('\t\t+')
                location_sample = location_transcript.resized( 0-sample_range       , 0-len(seq_transcript) )  # (shift to transcript_end, shift to sampling range)

            location_gene   = location_transcript
            seq_gene        = str(ensembl_genome.getRegion( location_transcript ).Seq)
            location_sample = location_sample
            seq_sample      = str(ensembl_genome.getRegion( location_sample     ).Seq)

        elif sample_direction=='downstream':
            transcript          = gene.getLongestCdsTranscript()
            seq_transcript      = str(transcript.Seq)
            seq_utr             = str(transcript.Utr3)
            location_transcript = transcript.Location

        # TAKE ANNOTAED UTR? If so is it AVAILABLE?
            # YES UTR: 
            if take_annotated_utr and len(seq_utr)>3:
                if test_it:
                    print('\t\tUTR annotated')
                annotated_sample_flag = True  # if UTR is annotated, no need to check for overlapping features
                utr_query      = re.compile(seq_utr)
                utr_match_all  = [m.start() for m in utr_query.finditer(seq_transcript)]
            # UTR matches to Transcript?
                # YES:
                if utr_match_all:
                    if test_it:
                        print('\t\t\tUTR Matched')
                    if not len(utr_match_all) == 1:
                        if test_it:
                            print('\t\t\t\tMultiple times')
                        sample_failures.append((geneId,'UTR matches multiple times to transcript...'))
                    utr_match_start = utr_match_all[-1] # Takes last matching instance as the start of match
                # NO:  Rare intronic transcript detected
                else:
                    if test_it:
                        print('\t\t\tUTR Not Matched')
                    #sample_failures.append((geneId,'UTR3 failed to match back to the transcript...'))

                    # CHOPPER MATCHER:  rare intronic exon prevents the annotated Utr from matching to it, so we chop the UTR from 5' to 3' until we can match back to seq_transcript
                    #
                    # seq_transcript:         -------TTG~~~TAG 
                    #
                    # annotated utr:                 TTGTAG     <= this cannot match to seq_transcript
                    #
                    # chopped utr:                   TTGTA >8   ...chopping away...
                    #                                TTGT >8  
                    #                                TTG >8     <= this can match batck!
                    #
                    utr_match_length = len(seq_utr)
                    while seq_utr[:utr_match_length] not in seq_transcript:
                        utr_match_length -= 1
                    seq_utr_chopped = seq_utr[:utr_match_length]
                    if len(seq_utr_chopped)<20:
                        sample_failures.append((geneId,'UTR is short yet does not match to transcript...'))
                    utr_query       = re.compile(seq_utr_chopped)
                    utr_match_all   = [m.start() for m in utr_query.finditer(seq_transcript)]
                    utr_match_start = utr_match_all[-1]

                if location_transcript.Strand == -1:
                    if test_it:
                        print('\t\t-')
                    location_utr                = location_transcript.resized( 0                    , 0-utr_match_start )  # (shift to utr_start, shift to utr_end)
                    location_transcript_noUtr   = location_transcript.resized( 0+len(seq_utr)       , 0 )                  # truncate the utr away
                else:
                    if test_it:
                        print('\t\t+')
                    location_utr                = location_transcript.resized( 0+utr_match_start    , 0 )                  # (shift to utr_start, shift to utr_end)
                    location_transcript_noUtr   = location_transcript.resized( 0                    , 0-len(seq_utr) )     # truncate the utr away
                seq_transcript_noUtr            = str(ensembl_genome.getRegion( location_transcript_noUtr).Seq)
            # NO UTR:
            else:
                if test_it:
                    print('\t\tUTR Synthetic')
                annotated_sample_flag = False
                if location_transcript.Strand == -1:
                    if test_it:
                        print('\t\t-')
                    location_utr                = location_transcript.resized( 0-sample_range       , 0-len(seq_transcript) )  # (shift to transcript_end, shift to sampling range)
                    location_transcript_noUtr   = location_transcript
                else:
                    if test_it:
                        print('\t\t+')
                    location_utr                = location_transcript.resized( 0+len(seq_transcript), 0+sample_range )  # (shift to transcript_end, shift to sampling range)
                    location_transcript_noUtr   = location_transcript
                seq_transcript_noUtr            = str(ensembl_genome.getRegion( location_transcript_noUtr ).Seq )
                seq_utr                         = str(ensembl_genome.getRegion( location_utr ).Seq )
        # MODE 1: If annotated UTR then keep it
            location_sample = location_utr 
            seq_sample      = seq_utr

            location_gene   = location_transcript_noUtr
            seq_gene        = seq_transcript_noUtr
        # MODE 2: Maintain Constant Sampling Length // Rob & Mike dont like it
            # location_sample = location_utr.resized(0,200-len(location_utr)) # e.g. if annotaed utr = 100bp, we add +100 = 200bp // or if utr = 300bp, we subtract 100bp = 200bp
            # seq_sample          = str(ensembl_genome.getRegion(location_sample).Seq)
            # gene_exons          = sorted(list(itertools.chain(*gene_exons)),key=attrgetter('End'))  # flattens from: 2d->1d
            # gene_limiting_exon  = gene_exons[-1]    # exon closest to 3' end of sample
            # location_sample     = gene_limiting_exon.resized(+len(gene_limiting_exon)+1,+sample_range+1) 

            sample_type = None
            if annotated_sample_flag:       # Labels the utr according to whether it was already annotated or synthesized by us
                sample_type = 'Annotated'
            else:
                sample_type = 'Synthetic'

    # INITIATE STORAGE:         ANDY: concatenate just the limiting exon sequence to the sampleseq?
        sample_seqs[geneId] = {'untruncated':str(),'truncated':str(),'gene_seq':str(seq_gene),'sample_location':str(),'gene_location':':'.join(str(location_gene).split(':')[2:]),'sample_type':sample_type}

    # FEATURES OVERLAP WITH SAMPLE?
        overlaps_flag = None 
        # NO    - since annotated samples means we assume no overlaps
        if annotated_sample_flag: 
            overlaps_flag = False
            location_sample_truncated = location_sample
        # MAYBE - it's not annotated, so let's check if our synthetic sample has overlaps...
        else:
            overlap_features = [feature for feature in ensembl_genome.getFeatures(region=location_sample,feature_types='gene') if not feature.StableId==geneId] # do any features overlap with our sample
            # YES but - so we've found overlapping genes, but specifically do any of their Exons overlap..?
            if overlap_features:        # If there are any overlaps we will need to enter truncation mode --> b) OVERLAPS
                overlaps_flag = True
                overlap_exons       = [[[exon for exon in transcript.Exons] for transcript in transcripts] for transcripts in [feature.Transcripts for feature in overlap_features]] # exons Of Transcripts Of Features
                overlap_exons       = list(itertools.chain(*list(itertools.chain(*overlap_exons)))) # flattens from: 3d->2d->1d
                overlap_locations   = [exon.Location for exon in overlap_exons]
                #================================================================================================================================================
                # PRUNE FALSE OVERLAPS:    remove exons that are so far down/upstream of sample they are beyond sample overlap // WARN: can be redundant?
                #================================================================================================================================================
                #
                #   overlapping exons:      5'----3'    5'----3'         5'----3'       5'----3' 
                #                              |||        |||||            xxxx           xxxx   		notice we dont want the xxx's to truncate away our Sample
                #   sample:                    5'-------S---------3'+5'-----------G------------------// S = sample, G = gene of sample, | = parts of sample to truncate
                #
                overlap_locations = [i for i in overlap_locations if (i.Start < location_sample.End) and (i.End > location_sample.Start)]   # remove "too-far-downstream" features

                                                                                                                                            # TODO put here? location_sample_truncated = location_sample
                # NO - the gene overlaps but not with any Exons of those genes!
                if overlap_locations==[]:   # If there are no viable overlaps (i.e. all overlaps were "beyond") then we keep the entire sample --> NO OVERLAPS
                    overlaps_flag = False
                    location_sample_truncated = location_sample                                                                         # TODO do we need an else here..?
            # NO - no genes overlap!
            else:
                overlaps_flag = False       # Else we keep the entire sample --> a) NO OVERLAPS
                location_sample_truncated = location_sample

    # SAMPLE TRUNCATION PROCEDURES:
        # NO OVERLAPS?  =>   Store whole sample, untouched...
        if not overlaps_flag:
            if test_it:
                print('\t\tRelevant Overlaps:  NO')
            sample_seqs[geneId]['untruncated']      = sample_seqs[geneId]['truncated'] = seq_sample
            sample_seqs[geneId]['sample_location']  = ':'.join(str(location_sample).split(':')[2:]) # genome coordinates summary
        # OVERLAPS?     =>   Truncate sample accordingly...
        else:
            if test_it:
                print('\t\tRelevant Overlaps:  YES')
            #===============================================================================================================
            # TRUNCATION STEP 1:  find the overlapping exon nearest to the sample-gene boundary, i.e. the "limiting exon"
            #===============================================================================================================
            #v--- a)/b) refers to the two main classes of routine for overlap truncation depending on whether ---v
            #   overlapping exons:                            5'---3' 5'---3' 5'---3'
            #                                                >>>|||<<<         |||
            #b) sample:                       //---G---3'+5'---S---3'          |||                      <=downstream +1 or upstream -1 sample
            # OR                                                            >>>|||<<<            
            #a) sample:                                                      5'---S---3'+5'---G---//    <=downstream -1 or upstream +1 sample
            #
            #   ...so after sorting we grab the limiting exon...
            #
            #b) limiting exon:                                5'---3' 
            # OR                                                                         ...these will then feed into STEP 2
            #a) limiting exon:                                                5'---3'
            #
            #===============================================================================================================
            # TRUNCATION STEP 2:  using the 'limiting exon' truncate the overlapping sample so it no longer ovelerlaps
            #===============================================================================================================
            #   limiting exon:                        5'------------3'      
            #                                           |||||   |||||  
            #a) sample:                                 |||||   5'---S---3'+5'---G---//                 <= downstream -1 or upstream +1 sample
            # OR                                        |||||
            #b) sample:                //---G---3'+5'---S---3'                                          <= downstream +1 or upstream -1 sample                                  
            #
            #   ...so after truncation we have...
            #   
            #a) sample:                              ---S....                   ... = truncated
            # OR
            #b) sample:                                          ....S---
            #

        # Which ORIENTATION do we truncate from: 
            # a) 3'-to-5' direction:                           <=i.e. gene is on the right hand side of the sample
            sample_completeOverlap_flag = None
            if (sample_direction == 'upstream' and location_sample.Strand == +1) or (sample_direction == 'downstream' and location_sample.Strand == -1):

                overlap_locations = sorted(overlap_locations,key=attrgetter('End'))   # Sort by 'End' <= 3' of feature vs. 5' of sample
                overlap_limit     = overlap_locations[-1]                             # limiting overlapping feature beyond which there are no overlaps w/i sample

            # Overlap SPANS entire sample?
                # YES:  truncate all
                if overlap_limit.End > location_sample.End:
                    if test_it:
                        print('\t\t\t\t\tcomplete overlap!')
                    sample_completeOverlap_flag = True
                    location_sample_truncated       = location_sample.copy()
                    location_sample_truncated.Start = location_sample_truncated.End
                # NO:   leave the remainder
                else:
                    if test_it:
                        print('\t\t\t\t\tpartial overlap!')
                    location_sample_truncated = location_sample.resized( 0+(overlap_limit.End-location_sample.Start) , 0 ) # Truncate the sampled upstream sequence (promoter to prevent overlap)
            # b) 5'-to-3' direction:                           <=i.e. gene is on the left hand side of the sample
            elif (sample_direction == 'downstream' and location_sample.Strand == +1) or (sample_direction == 'upstream' and location_sample.Strand == -1):

                overlap_locations   = sorted(overlap_locations,key=attrgetter('Start')) # 'End' <= 3' of feature vs. 5' of sample
                overlap_limit       = overlap_locations[0]

                if overlap_limit.Start < location_sample.Start:
                    if test_it:
                        print('\t\t\t\t\tcomplete overlap!')
                    sample_completeOverlap_flag = True
                    location_sample_truncated       = location_sample.copy() 
                    location_sample_truncated.End   = location_sample_truncated.Start
                else:
                    if test_it:
                        print('\t\t\t\t\tpartial overlap!')
                    location_sample_truncated = location_sample.resized( 0 , 0-(location_sample.End-overlap_limit.Start) ) # Truncate the sampled upstream sequence (promoter to prevent overlap)

            # if location_sample.Strand == -1:
            #     return overlap_features, overlap_limit, overlap_locations, location_sample, location_transcript, location_transcript_noUtr, ensembl_genome

        # HOW MUCH OVERLAP?
            # COMPLETE
            if sample_completeOverlap_flag:
                sample_seqs[geneId]['untruncated']      = str(ensembl_genome.getRegion(location_sample).Seq)            # untruncated asmple sequence
                sample_seqs[geneId]['truncated']        = ''                                                            # No sample sequence due to complete overlap
                sample_seqs[geneId]['sample_location']  = ':'.join(str(location_sample_truncated).split(':')[2:])       # genome coordinates summary
            # PARTIAL
            else:
                sample_seqs[geneId]['untruncated']      = str(ensembl_genome.getRegion(location_sample).Seq)        
                sample_seqs[geneId]['truncated']        = str(ensembl_genome.getRegion(location_sample_truncated).Seq)  # truncated sample sequence
                sample_seqs[geneId]['sample_location']  = ':'.join(str(location_sample_truncated).split(':')[2:])      

    # MASKING PROCEDURE?
        # YES
        if not masking=='none':

            if overlaps_flag:
                if sample_completeOverlap_flag:
                    continue
                else:
                    location_sample = location_sample_truncated
            else:
                location_sample = location_sample

            sample_seq          = list(sample_seqs[geneId]['truncated'])
            repeats             = ensembl_genome.getFeatures(region=location_sample,feature_types='repeat')
            repeats_regions     = [[i.Location.makeRelativeTo(location_sample).Start,i.Location.makeRelativeTo(location_sample).End] for i in repeats]

            # MASK against each repeat:
            for region in repeats_regions:
                # Are repeat.Location's start or end out of bounds? --> restrict them to sample_locations limits if so
                if region[0]<0:
                    region[0]=0
                if region[1]>len(sample_seq):     # TODO: is len() the actual matcher? 
                    region[1]=len(sample_seq)
                if masking=='hard':
                    sample_seq[region[0]:region[1]+1]='N'*(len(sample_seq[region[0]:region[1]+1])) # PyCogent is not happy with X
                elif masking=='soft':
                    sample_seq[region[0]:region[1]+1]=[letter.lower() for letter in sample_seq[region[0]:region[1]+1]]
            sample_seqs[geneId]['truncated']=''.join(sample_seq)

            #print(sample_seq)

        # NO
        # ...then move on...

    return sample_seqs, sample_direction, sample_failures  # yields the dictionary      # ANDY_02/03: 

#______________________________________________#
# write sequences into either a pickle or fasta
#______________________________________________#
def sampleSequences_write(ensembl_genome, sample_read, fasta_it=True, pickle_it=True, include_genes=False):
    """ 
    Notes:
        Stores into files (either .fasta or .p) the sampled sequence data that was generated into a dict by sampleSequences_read()

    Args:
        ensembl_genome  = setupGenome( ... )                    # ensembl genome object: http://pycogent.org/examples/query_ensembl.html
        sample_read     = sampleSequences_read( ... )       # tuple: ( dict of sampled seqs, sample_direction e.g. 'upstream')
        fasta_it        = True                                                  # bool: generate a fasta file for all seqs? 
        pickle_it       = False                                              # bool: store a python dict with all the sequences as file?
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
            total_seq       = gene_seq+sample_seq 
            gene_length     = len(gene_seq)
            sample_length   = len(sample_seq)
            header = geneId+'\tGeneCoord:'+sample_data[geneId]['gene_location']+'\tUtrCoord:'+sample_data[geneId]['sample_location']+'\tGeneLength:'+str(gene_length)+'\tUtrLength:'+str(sample_length)+'\tUtrType:'+str(sample_data[geneId]['sample_type'])
        else:
            total_seq       = sample_data[geneId]['truncated']
            total_length    = len(total_seq)
            header = geneId+'\tUtrCoord:'+sample_data[geneId]['sample_location']+'\tUtrLength:'+str(total_length)+'\tUtrType:'+str(sample_data[geneId]['sample_type'])

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
        samplelength_distribution = [int(i[4].replace('UtrLength:','')) for i in headers]
    else:
        headers = [i.split('\t') for i in samples.keys()]
        samplelength_distribution = [int(i[2].replace('UtrLength:','')) for i in headers]
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
            genome              = setupGenome(              species, db_host='localhost', db_user='vbuser', db_pass='Savvas', db_release=73 )
            samples_read        = sampleSequences_read (    genome, sample_direction=sample_direction, sample_range=500, take_annotated_utr=True, masking = 'hard', sample_data={}, sample_seqs={}, test_it=False)
            samples_write   = sampleSequences_write(        genome, samples_read, fasta_it=True, pickle_it=True, include_genes=False)
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

""" NOTICE:

    - 05/02: sample_range 		= 2000
    - 06/02: take_annotated_utr = False
    - 10/02: masking            = 'hard'

"""

# sample_directions = ['upstream']
# sampleAllSpecies(species_list,sample_directions)

# sample_directions = ['downstream']
# sampleAllSpecies(species_list,sample_directions)

#-------------------------------------------------------------------------------------------------------------------------------------
# @TESTING
#-------------------------------------------------------------------------------------------------------------------------------------

# genome          = setupGenome(              'Anopheles gambiae', db_host='localhost', db_user='vbuser', db_pass='Savvas', db_release=73 )
# samples_read    = sampleSequences_read (    genome, sample_direction='upstream', sample_range=2000, take_annotated_utr=False, masking = 'none', sample_data={}, sample_seqs={}, test_it=True)
# samples_write   = sampleSequences_write(    genome, samples_read, fasta_it=True, pickle_it=True, include_genes=False)

# #overlap_features, overlap_limit, overlap_locations, location_sample, location_transcript, location_transcript_noUtr

