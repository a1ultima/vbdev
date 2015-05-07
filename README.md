vbdev
=====

VB Developer Backup

#####################################################
# README v3 2014-11-18 v2 2014-11-14  v1 2014-11-05 #
#####################################################

BEFORE READING

    This is a README to help understand the scripts and data from ./, for more information on the data that is related to the scripts see /home/ab108/0VB/2kb/data/README. 

    Often you can find much more specific descriptions within the headers of each script, so called "MINI README". Also, you will sometimes see more fine grained info within READMEs in sub-directories. 

    If there are special characters within: (i) READMEs, (ii) MINI READMEs (header section of scripts) or (iii) Comments within scripts that do not make sense to you see the "FORMATTING STANDARDS KEY:" section of /home/ab108/0VB/2kb/README. If there are words that do no make sense to you but have a leading @ symbol then see the "GLOSSARY:" section of /home/ab108/0VB/2kb/README.

    Much of the statistics that were shared with VB is available: https://docs.google.com/spreadsheets/d/1TVAtKTi6Rg6u_RJrv00WGd3U8SeMjl3xAuaakrPV0lk/edit#gid=320127454, and work logs for this project is available: https://docs.google.com/document/d/1WMyT3leCDB02UYzmhg6M3uIlDLf8l0Vu7y-PKHvHFZE/edit# feel free to ask for permission to view them.


DESCRIPTION: 

    Scripts here are responsible for (i) preparing "promoteromes" ** of 21 newly sequenced/assembled species ***, (ii) discovering motifs using @dreme., (iii) Merging motifs discovered in separate species into "regulatory motif families", (iv) Determining where these motifs are distributed across the promoterome.

    ** 2kb DNA fragments that are upstream of genes of the genomes of species specified by ./scripts/species_list.txt (for more info see comments within: ./updownstream.py).

    *** you can see the list of species whose genomes we've discovered motifs for at: home/ab108/0VB/2kb/scripts/species_list.txt.



PIPELINE MAP: 

@Species included in the pipeline can be seen as the non-commented species names specified by ./scripts/species_list.txt

Phase 1: Discover Motifs 

    all of PHASE 1 is coordinated by: /home/ab108/0VB/2kb/scripts/pipeline_meme.py

    p1.1. python /home/ab108/0VB/2kb/scripts/updownstream.py: 
        
        Generate a set of "promoteromes", one per 21 newly sequenced & assembled species avialable in the VectorBase EnsEMBL database (see: /home/ab108/0VB/2kb/scripts/species_list.txt).

        data in:

            mysql -hlocalhost -uvbuser -pSavvas # EnsEMBL genome databse hosted on vb-dev.bio.ic.ac.uk

            /home/ab108/0VB/2kb/scripts/species_list.txt # config file to specify which species to include in the pipeline

        data out: 

            /home/ab108/0VB/2kb/data/sample_seqs/fasta/* # where * is <@species>_upstream.fasta see python /home/ab108/0VB/2kb/scripts/updownstream.py for more info.


    p1.2. python /home/ab108/0VB/2kb/scripts/meme_dataPrepper.py 
        
        Prepare @dreme input data (promoteromes) for each species subject to dreme requirements. Mainly: performs dustmasking, removes redundant sequences, then formats promoteromes so that each line has only 100bp. 

        data in:

            /home/ab108/0VB/2kb/scripts/species_list.txt
            /home/ab108/0VB/2kb/data/meme_data/in/

        data out:

            /home/ab108/0VB/2kb/data/meme_data/in/<@species>_upstream_memeready_all.fasta
            /home/ab108/0VB/2kb/data/meme_data/in/<@species>_upstream_memeready_all_simpleMasked.fasta

    p1.3. python /home/ab108/0VB/2kb/scripts/dreme_randSampleFasta.py: 

        @todo (see MINI README within script)

        data in: 

            /home/ab108/0VB/2kb/scripts/species_list.txt
            /home/ab108/0VB/2kb/data/meme_data/in/<@species>_upstream_memeready_all_simpleMasked.fasta

        data out:

            /home/ab108/0VB/2kb/data/meme_data/in/random_dreme/*

    p1.4.a. Single Core (vb-dev): python /home/ab108/0VB/2kb/scripts/dreme_params.py  OR 
    p1.4.b. Multi Core (HPC): bash /home/ab108/0VB/2kb/scripts/hpc/dreme_submit_all_species.sh

        Run @dreme for all species. Note: Culex quinquefasciatus, and Aedes aegypti take >9 days on vb-dev.bio.ic.ac.uk. One core dedicated to each. The others finish in ~3 days and can be done on HPC (one core per species).

        ETA: ~3 days per species. 9+ days for Aedes aegypti or Culex quinquefasciatus.

        data in:

            /home/ab108/0VB/2kb/data/meme_data/in/** 

                ** <@species>_upstream_memeready_all_simpleMasekd.fasta e.g. anopheles_gambiae_upstream_memeready_all_simpleMasked.fasta

        data out:

            Single Core: /home/ab108/0VB/2kb/data/meme_data/out/dreme_100bp/evalue_0.05/**  OR
            Multi Core (HPC): /home/ab108/0VB/2kb/data/meme_data/out/dreme_100bp/sampled_all_hpc/evalue_0.05/** 

                ** <species>_100bp/*** e.g. anopheles_gambiae_100bp

                    *** dreme.fasta, dreme.txt, dreme.html, dreme.xml, <motifname>.png

    p1.5. { Manual Steps, see: /home/ab108/0VB/2kb/scripts/manual/* } 

        Manual steps just involve: (i) placing files needed by STAMP in the right place and (ii) shipping @dreme discovered motifs to be located in /home/ab108/0VB/2kb/data/meme_data/out/@dreme_100bp/sampled_all_hpc/evalue_e0.05/. If they are not present in that location, move them there. The data should appear like so: drosophila_melanogaster_100bp/@dreme.txt anopheles_gambiae_100bp/@dreme.txt, etc. If you are really stuck, then have a look at ./manual/stamp.sh, which contains several manual command-line snippets used at the time. If that doesn't help maybe have a look at ./manual/stamp_prepper.sh. If that still doesnt help, email me: andrew.i.brockman@gmail.com!

        data in:

            N/A

        data out:

            /home/ab108/0VB/2kb/data/stamp_data/in/common/

    p1.6. bash /home/ab108/0VB/2kb/scripts/dreme_to_stamp.sh 

        After dreme has discovered motifs for each of the @species this script will prepare**** the data ready for @STAMP***** to generate a motif tree*** then run @STAMP. 

            Please take note of: ctrl+f: "ECUT=", this specifies the e-value cut-off that dreme used to discover motifs. The standard we went with was that e-value = 0.05
            . That is, ECUT=0.05 

            *** Preparation involves (i) converting motif @PWMs discovered separately for each @species from .meme format into STAMP-compatible TRANSFAC format, (ii) uniting the PWMs in separate directories/files into a single file: /home/ab108/0VB/2kb/data/meme_data/out/dreme_100bp/sampled_all_hpc/evalue_0.05/all_100bp/dreme.fasta. 

            **** The motif tree consists of motifs organised hierarchically by similarity of aligned @PWMs. The tree allows us to merge motifs discovered separately for each @species into "motif families". The tree can then be collapsed at parent nodes which will merge together motifs from separate species that are functionally alike (orthologues / paralogues). 

            ***** see: http://www.benoslab.pitt.edu/stamp/. As for the arguments we chose for STAMP, -cc SSD -align SWU -out were chosen based on Andy and Bob's eyeballing of how the trees looked after testing out all options. 

        data in: 

            /home/ab108/0VB/2kb/scripts/species_list.txt
            /home/ab108/0VB/2kb/data/stamp_data/in/common/jaspar.motifs # stamp:required:file:@jaspar database of motifs in TRANSFAC format which is targeted by the input motifs to STAMP for matching

        data out: 

            /home/ab108/0VB/2kb/data/stamp_data/in/dreme_100bp_e$ECUT/dreme.fasta" 
            /home/ab108/0VB/2kb/data/stamp_data/in/common/ScoreDists/JaspRand_SSD_SWU.scores"
            /home/ab108/0VB/2kb/data/stamp_data/out/dreme_100bp_e$ECUT/SWU_SSD/out"

Phase 2: Analyse Motifs (coordinated by: /home/ab108/0VB/2kb/scripts/stamp_to_collapsed/...):

    p2.1. python /home/ab108/0VB/2kb/scripts/stamp_to_collapsed/collapseMotifTree.py

        generate a tree of motif similarity and collapse it at various depths, each collapsed node is a set of merged motifs, a motif cluster 

        WARNING: see commenting available in the script, try ctrl+f: "@WARN" within the script.

        ETA: 36 hrs 
    
    p2.2. python /home/ab108/0VB/2kb/scripts/stamp_to_collapsed/motifStatistics.py

        <describe it @todo>

        data in: 

            /home/ab108/0VB/2kb/data/stamp_data/out/* # everything except the directory:/progressively_collapsed_motifs/*

        data out: 



        Generates an essential datastructure to downstream analyses /home/ab108/0VB/2kb/data/stamp_data/out/@dreme_100bp_e0.05/SWU_SSD/cluster_to_stats.p, which is essentially a giant python dict storing various statistics both at the per @dreme motif level as well as the "motif @cluster" level. 

        WARNING: see commenting available in the script, try ctrl+f: "@WARN" within the script.

        ETA: 10 hrs (?)
 
    p2.3. python /home/ab108/0VB/2kb/scripts/stamp_to_collapsed/collapseMotifTree_progressiveMode.py : blacklist_then_summaryStats_S_and_H_combos( S_from = 1, S_to = 21, S_step = 1, H_from = 1, H_to = 30, H_step = 1, e = 0.05)

        <describe it @todo>

        For more detailed information see the mini README within the script.

        ETA: 36 hrs
 
    p2.4. cd /home/ab108/0VB/2kb/scripts/stamp_to_collapsed/progressiveMode_to_tomtom/ python format_db_cisbp_to_meme.py --> bash combos_to_tomtom_master.sh
        
        comparing our motifs to target databases such as: CIS-BP, JASPAR CORE, JASPAR INSECTS and various Drosophila databases: OnTheFly_2014_Drosophila.meme, dmmpmm2009.meme, flyreg.v2.meme, idmmpmm2009.meme.
        
        p2.4.1 bash vennDiagram_master.sh
        
            Firstly it converts CIS-BP motifs to .meme format. This then allows us to make comparisons of our motifs to CIS-BP motifs using TOMTOM. Then it unites all unique CIS-BP motifs that have a match to a transcription factor listed in /home/ab108/software/meme/db/CIS-BP/downloaded_format/*/TF_Information.txt, unites them into a single .meme file and runs TOMTOM with this .meme file as a target and query:@AGCC-final-motif-set

        p2.4.2

        @dropped:

            combos_to_tomtom_master.sh

                Runs tomtom using query: "dipteran" or query: "all" vs. target databases: "JASPAR CORE 2014, JASPAR INSECTS 2014, CIS-BP, etc". Then calls upon another script to calculate statistics of the fraction of our motifs match to target databses, in files such as: /home/ab108/0VB/2kb/data/stamp_data/out/@dreme_100bp_e0.05/SWU_SSD/progressively_collapsed_motifs/nspecies_3.0_entropy_5.0/tomtom/match_counts.txt

            WARNING: there is a parameter that controls whether it runs on single or parallel. If you want it to run on parallel you must manually edit the scripts in the three scripts in /home/ab108/0VB/2kb/scripts/stamp_to_collapsed/progressiveMode_to_tomtom/parallel/. So that they iterate through the intented parameterisations selected in Phase 2: 3.'s script. If in doubt just set it to "single".

        ETA: 48 hrs

    p2.5. bash /home/ab108/0VB/2kb/scripts/collapsed_to_fimo/tomtom_to_fimo_master.sh 

        Runs FIMO to find occurrences of motifs in the promoteromes of all @species (see: /home/ab108/0VB/2kb/scripts/collapsed_to_fimo/README ).

        ETA: 36 hr (?) per @species.




GLOSSARY (and info  @todo):

@CLUSTER: 

    a motif cluster is a group of motifs that are merged together. Merging means the combining of each constituent motif's Position Weight Matrix (@PWM). In terms of the motif tree, a motif cluster is the parent node that joins children nodes, where the children nodes are individual motifs. 

@MOTIF: 

    a motif that was originally discovered by @dreme when fed a @species's @promoterome as input data. 

@SPECIES

    one of the 21 newly sequenced/assembled insect species we are interested in. 

    Not all species of the 21 available were used. For a list of species included in this motif discovery pipeline see the non-commented** lines in ./scripts/species_list.txt

    Around 2014-11-xx, @Andy, @Bob and @Giannis decided it might be best to drop Aedes aegypti and Culex quinquefasciatus.

    **commented: lines of files prefixed with #

@PROMOTEROME: 

    is the sequence data generated by 2kb/scripts/upDownStream.py, which are DNA sequences upstream of every gene per non-commented EnsEMBL genome species in /home/ab108/0VB/2kb/scripts/species_list.txt. Upstream of every gene, DNA sequences of length N is taken if available, otherwise 2kb is sampled ad hoc, where N is the length of 5'UTR already annotated by EnsEMBL.

@PWM: 

    Position Weight Matrix, a.k.a Family Binding Profile (FBP), a.k.a PSPM, a.k.a. PSSM is the form of data that motifs are encoded in. Each @ROW represents a position along a motif sequence. The four decimals in each row are relative frequencies of A C G T letters in that order (e.g. leftermost decimal in a row is A and rightermost is T, note: in alphabetical order).

    See: Wiki: "PWM". 

    These come in various formats: TRANSFAC, .meme, etc. For examples of .meme 

            see: /home/ab108/0VB/2kb/data/stamp_data/out/@dreme_100bp_e0.05/SWU_SSD/progressively_collapsed_motifs/distance_cut_0.1/nspecies_3.0_entropy_5.0/all_meme_format.txt, for examples of TRANSFAC 

            see: /home/ab108/software/meme/db/CIS-BP/transfac_format/*

@DREME: 

    from http://meme.nbcr.net/meme/doc/dreme.html: "DREME (Discriminative Regular Expression Motif Elicitation) finds relatively short motifs (up to 8 bases) fast, and can perform discriminative motif discovery if given a negative set, consisting of sequences unlikely to contain a motif of interest that is however likely to be found in the main ("positive") sequence set. If you do not provide a negative set the program shuffles the positive set to provide a background (in the role of the negative set). The input to DREME is one or two sets of DNA sequences. The program uses a Fisher Exact Test to determine significance of each motif found in the postive set as compared with its representation in the negative set, using a significance threshold that may be set on the command line. DREME achieves its high speed by restricting its search to regular expressions based on the IUPAC alphabet representing bases and ambiguous characters, and by using a heuristic estimate of generalised motifs' statistical significance." 

    EXAMPLE OUTPUT see: /home/ab108/0VB/2kb/data/meme_data/out/dreme_100bp/sampled_all_hpc/evalue_0.05/anopheles_gambiae_100bp/dreme.txt (note other outputs include *.png motif logos, dreme.html and dreme.xml)

    EXAMPLE INPUT see: /home/ab108/0VB/2kb/data/meme_data/in/random_dreme
/anopheles_gambiae_upstream_dremeready_all_simpleMasked_random.fasta. 

@FIMO: 
    
    from http://meme.nbcr.net/meme/doc/fimo.html: "The name FIMO stands for "Find Individual Motif Occurences." The program searches a database of DNA or protein sequences for occurrences of known motifs, treating each motif independently. The program uses a dynamic programming algorithm to convert log-odds scores (in bits) into p-values, assuming a zero-order background model. By default the program reports all motif occurrences with a p-value less than 1e-4. The threshold can be set using the --thresh option. The p-values for each motif occurence are converted to q-values following the method of Benjamini and Hochberg ("q-value" is defined as the minimal false discovery rate at which a given motif occurrence is deemed significant). The --qv-thresh option directs the program to use q-values rather than p-values for the threshold. If a motif has the strand feature set to +/- (rather than +), then fimo will search both strands for occurrences. The parameter --max-stored-scores sets the maximum number of motif occurrences that will be retained in memory. It defaults to 100,000. If the number of matches found reaches the maximum value allowed, FIMO will discard 50% of the least significant matches, and new matches falling below the significance level of the retained matches will also be discarded. FIMO can make use of position specific priors (PSP) to improve its identification of true motif occurrences. To take advantage of PSP in FIMO you use must provide two command line options. The --psp option is used to set the name of a MEME PSP file, and the --prior-dist option is used to set the name of a file containing the binned distribution of priors." 

    EXAMPLE input see: 
        input1 (@promoterome): /home/ab108/0VB/2kb/data/fimo/in/fasta/anopheles_gambiae_upstream.fasta, 
        input2 (@motif file): /home/ab108/0VB/2kb/data/stamp_data/out/dreme_100bp_e0.05/SWU_SSD/progressively_collapsed_motifs/distance_cut_0.35/nspecies_1.0_entropy_30.0
/all_meme_format.txt

@DROPPED:

    a script or data that is no longer necessary for the pipeline. At the time these were created during the exploratory periods of the projects.

@STAMP:

    Tool for aligning @motifs for comparison, it also generates newick @trees of both the entire set all off @species' motifs generated by @dreme as well as 
    
    for more info:

        type "stamp" in shell.

        see: http://www.benoslab.pitt.edu/stamp/. 

    arguments/options we chose in pipeline:

        As for the arguments we chose for STAMP, -cc SSD -align SWU -out were chosen based on Andy and Bob's eyeballing of how the trees looked after testing out all options. 

@TREE:

    The motif tree consists of motifs organised hierarchically according to similarity in aligned motif Position Weight Matrices (PWMs). The tree allows us to unite and comapre motifs discovered individually for each @species into a single data structure. The tree can then be collapsed at parent nodes which will merge together motifs from separate species that are functionally alike (orthologues / paralogues). 

    Format is always in Newick:type: (:0.1,:0.2,(:0.3,:0.4):0.5):0.0;       all have a distance to parent

    for more info:

        see:wiki: http://en.wikipedia.org/wiki/Newick_format


@AGCC-final-motif-set:

    135 motifs decided to be our final set of motifs. These were chosen on the basis that their distribution of entropies was the most similar to CIS-BP's Anopheles gambiae, Aedes agypti and Culex quinquefasciatus distribution of entropies.

    PWMs of the final set can be found:

        /home/ab108/0VB/2kb/data/stamp_data/out/dreme_100bp_e0.05/SWU_SSD/progressively_collapsed_motifs/distance_cut_0.35/nspecies_1.0_entropy_30.0/all_meme_format.txt

    After making comparisons between CIS-BP and AGCC-final-motif-set, we took the AGCC-specific motifs + AGCC motifs with a match to CIS-BP + CIS-BP specific motifs to generate a AGCC and CIS-BP set:

        /home/ab108/0VB/2kb/data/stamp_data/out/dreme_100bp_e0.05/SWU_SSD/progressively_collapsed_motifs/distance_cut_0.35/nspecies_1.0_entropy_30.0/cisbp_agcc_final_merged_set.meme




FORMATTING STANDARDS KEY:

Formatting standards apply to: READMEs, MINI READMEs (header section of a script), Comments (of a script) and to any of the following bodies of text (delimited by "standards") IF there is a standard already set above it.

<word>:

    Words surrounded by < > describe the general case of possible words that are valid within the < >. I.e. the literal word used in <word> is in some higher ontological level than possible words that are valid replacements of <word>. 

    For example if I want to refer to a file TF_Information.txt that is shared between several species I can refer to it's path as /home/ab108/software/meme/db/CIS-BP/downloaded_format/<species>/TF_Information.txt.

* (single): 
    
    same meaning as it does in most regular expression standards. The asterix is often there to be perceived by you as "general case", in less often cases it may also mean "and so on" (e.g. 1, 2, 3, *).

    E.G. if I am describing a type of data that is available for various directories whose names relate to input parameters of some program, I may refer to the path of such data as /<somepath>/*/<somedata>.txt. 

    It is for you to go and cd /<somepath>/ to see what files/sub-directories I was referring to when using *.

** (multiple: e.g. **, ***, ****):

    more information is available for a word or phrase that is suffixed with ** nearby and below. More information associated with a word is then prefixed with the same number of asterices. The available more information associated with the word suffixed by ** is only valid for the nearest matching **, i.e. there can be multiple occurrences of ** or *** or **** in this README but only the nearest pair of ** refer to one another.

@word:

    A word with a description available in the "GLOSSARY:" section. Words in the "GLOSSARY:" section are capitalised to help with ctrl+f efficiency.

@todo:

    A tag to indicate there is more for @Andy to do. 

    If you encounter this and you need more infomation or code to be developed in place of the @todo contact @Andy (ctrl+f:@Andy).

<word1>:<word2>:<word3>:*

    Words further to the left of : are higher up in an ontological hierarchy that makes sense somewhere in @Andy's brain.

@<first name of person>:

    Tags of names of people involved in the pipeline. And to help the reader to track contact details. E.g. @Bob, @Giannis, @Andy.

    @BOB:email:uncoolbob@gmail.com or r.maccallum@imperial.ac.uk
    @BOB:status:2014-11-28:line manager

    @GIANNIS:email:i.kirmitzoglou@imperial.ac.uk
    @GIANNIS:status:2014-11-28:research associate (?) 

    @ANDY:email:andrew.i.brockman@gmail.com or ab108@ic.ac.uk or a1ultima@gmail.com
    @ANDY:status:2014-11-28:research assistant

@param:

    only found in scripts, a tag to indicate that this line of code involves varaible assignment which can be changed as you would a standard input for a bash shell command. Often it has not been made into a proper standard input variable because I was prioritising time elsewhere.



