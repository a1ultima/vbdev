import updownStream 
import meme_dataPrepper
import meme_bgfileGen

speciesFile = './species_list.txt'
seqDataPath = '../data/sample_seqs/fasta/'

updownStream.sampleAllSpecies(  speciesFile         = speciesFile,
                                sample_directions   = ['upstream'],
                                dataPath_out        = '../data/sample_seqs/fasta/'
                                )
meme_dataPrepper.allSpecies(    speciesFile         = speciesFile, 
                                dataPath_in         = 'home/ab108/0VB/2kb/data/sample_seqs/fasta/', 
                                dataPath_out        = 'home/ab108/0VB/2kb/data/meme_data/in/',
                                LCR_masking         = 'simple'
                                )
meme_bgfileGen.allSpecies(      dataPath_in         = '0VB/2kb/data/meme_data/in/',
                                bfileGeneratorPath  = '0VB/2kb/scripts/', 
                                maskingChar         = 'n', 
                                order               = 3, 
                                speciesFile         = speciesFile
                                )
