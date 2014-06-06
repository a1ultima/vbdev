import updownStream 
import meme_dataPrepper
import meme_bgfileGen
import meme_randSampleFasta
import dreme_randSampleFasta
import dreme_params
import tomtom_params

speciesFile     = './species_list.txt'
dataPath_fasta  = '../data/sample_seqs/fasta/'
dataPath_meme   = '../data/meme_data/in/'

def configCheck(): # Recursively check if the user has properly configured the species_list.txt
      choice = raw_input("Did you make sure to configure ./species_list.txt ??? (y/n)")
      if choice=='y':
            print("Okay!")
      elif choice=='n':
            raise Exception("Thought so... here this'll help: vi ./species_list.txt")
      else: 
            configCheck()
configCheck() #^

updownStream.sampleAllSpecies(      speciesFile             = speciesFile,
                                    sample_directions       = ['upstream'],
                                    dataPath_out            = dataPath_fasta,
                                    sample_range            = 1000,  # 14-05-11: 2k bp -> 1k bp <= HPC takes > 72hr and at least this will reduce overlap pulley seqs
                                    take_annotated_utr      = False,
                                    masking                 = 'none' # 14-05-11: turns out this didnt affect very much afterall --| 14-05-10: make sure this shouldn't be set to "simple" because "simple" it may in some way have been an implicit arg...
                                    )

meme_dataPrepper.allSpecies(        speciesFile             = speciesFile, 
                                    dataPath_in             = dataPath_fasta,       
                                    dataPath_out            = dataPath_meme,
                                    LCR_masking             = 'simple'
                                    )

# meme_bgfileGen.allSpecies(          dataPath_in             = dataPath_meme,
#                                     bfileGeneratorPath      = '../scripts/', 
#                                     maskingChar             = 'n', 
#                                     order                   = 3, 
#                                     speciesFile             = speciesFile
#                                     )

# meme_randSampleFasta.allSpecies(    speciesFile             = speciesFile, 
#                                     dataPath_in             = dataPath_meme, 
#                                     dataPath_out            = '../data/meme_data/in/randomFasta/', 
#                                     n_replicates            = 5, 
#                                     n_seqs                  = 2000 
#                                     )

dreme_randSampleFasta.allSpecies(   speciesFile         = speciesFile,
                                    dataPath_in         = dataPath_meme,
                                    dataPath_out        = dataPath_meme+'random_dreme/',
                                    n_seqs              = 100000,
                                    len_seq             = 100,
                                    population          = True
                                    )

# dreme_params.allSpecies(            speciesFile         = speciesFile,
#                                     dataPath_in         = dataPath_meme+'random_dreme/',
#                                     dataPath_out        = '../data/meme_data/out/dreme_100bp/sampled_all/',
#                                     verbosity           = 2,
#                                     eValue              = 0.0001,
#                                     nMotifs             = 1000,
#                                     resultFormat        = 'png',
#                                     )

# tomtom_params.allSpecies(           speciesFile       = speciesFile,
#                                     dataPath_in       = '../data/meme_data/out/dreme_100bp/sampled_all_hpc/',
#                                     dataPath_out      = '../data/meme_data/out/tomtom_100bp/',
#                                     verbosity         = 1,
#                                     threshold         = 0.1
#                                     )
