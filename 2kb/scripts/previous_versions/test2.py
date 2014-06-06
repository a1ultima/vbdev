




meme_dusted = [     'meme', 
                    'anopheles_gambiae_upstream_memeready_all_dusted_100.fasta', 
                    '-dna', 
                    '-mod', 'oops', 
                    '-revcomp', 
                    '-w', '10', 
                    '-nmotifs', '5', 
                    '-maxsize', '50000000', 
                    '-p', '4', 
                    '-oc', '../out/ag'
                    ]




def memeScaler( argsList, argToVary_choice, argToVary_parameterSpace ):

    """ 
    Notes:
        Allows one to run MEME several times each with varying parameterisations

    Args:
        argToVary_choice = 7            # index to the arg you wish to vary, e.g. 7 => -w 10, where 10 will be varied

    """

    #n = len(argToVary_parameterSpace)

    outfile_prefix  = '../out/ag'   # TODO: can use part of the infile name as prefix
    outfile_names   = [outfile_name + argsList[7-1] + str(i) for i in argToVary_parameterSpace]

    for j,value in enumerate(argToVary_parameterSpace):

        argsList[-1] = outfile_names[j]     # alter the outfile name according to new parameterSace value

        for i,arg in enumerate(argsList):

            if i == argToVary_choice:

                argsList[i] = value # the argList's value is now changed according to parameterSpace

                start   = subprocess.Popen(['date'], stdout=subprocess.PIPE).communicate()[0]
                subprocess.call(argsList)    # Run meme with new argToVary now changed
                end     = subprocess.Popen(['date'], stdout=subprocess.PIPE).communicate()[0]
                tdelta  = datetime.strptime(end.split(' ')[3], '%H:%M:%S') - datetime.strptime(start.split(' ')[3], '%H:%M:%S')
                print(tdelta)