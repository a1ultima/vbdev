
"""

1. Read in the lines of "new" ~/popbio/data_andy/isatabs/andy_0411_MAP/old/Europe/a_collection.txt and "new" ~/0VB/isatab/malariaAtlas/data/isatab/Europe/a_collection.txt 

2. Sort the lines of "old" and "new" 

3. Check for equality

"""

#===============#
#   Imports     #
#===============#

import os 
import numpy as np 


#===============#
#   Functions   #
#===============#

def abs_to_rel( path ):
    """

    Description:

        Returns a relative path given the present working directory 

    Arguments:

        path         # string that is the absolute path 

    Usage:

        datapath = rel_to_abs('/home/ab108/0VB/isatab/malariaAtlas/scripts/test')
        fi = open( datapath )
        while True:
            line = fi.readline()
            if line == '':
                break
            print line
        fi.close()

    """
    import os 
    return os.path.relpath( path, os.getcwd() )


def diff(path1,path2):

    # Load old data (Europe a_collections.txt)

    old_path = abs_to_rel(path1)
    old_file = open(old_path)
    old_rows = []

    while True:
        line = old_file.readline()
        if line == '':
            break
        old_rows.append(line)
    old_file.close()

    # Load new (dates fixed) data (Europe a_collections.txt)

    new_path = abs_to_rel(path2)
    new_file = open(new_path)
    new_rows = []

    while True:
        line = new_file.readline()
        if line == '':
            break
        new_rows.append(line)
    new_file.close()


    new_diff = [row for row in sorted(new_rows) if row not in old_rows]
    old_diff = [row for row in sorted(old_rows) if row not in new_rows]

    for i in range(0,len(old_diff)):
        print '<'+old_diff[i]
    print '---'
    for i in range(0,len(new_diff)):    
        print '>'+new_diff[i]







#---------------------------------------------------------
# Scraps


# path1 = '/home/ab108/popbio/data_andy/isatabs/andy_0411_MAP/old/Europe/a_collection.txt'
# path2 = '/home/ab108/0VB/isatab/malariaAtlas/data/isatab/Europe/a_collection.txt'

# diff(path1,path2)


#---------------------------------------------------------



import sys,getopt


#@WARN: this script can be run from termnial, but NEEDS two arguments to be fed to it: (i) from terminal: --nspecies 3.0 --entropy 10.0, (ii) from Python: blacklist, summary, clade_groups, cluster_to_stats = blacklist_then_summaryStats_from_shell( e = 0.05, species_threshold = 3.0, entropy_threshold = 10.0 )

# Argument Handling
def main(argv):

    # USER MADE ERROR? 

    # ... if no
    try:

        opts, args = getopt.getopt(argv,"hf:s",["file1=","file2="])

    # ... if yes (e.g. mispelled args, non-existant args etc)
    except getopt.GetoptError:
        
        print 'diff.py --file1 <path to first file> --file2 <path to second file>' # <-- @@@

        sys.exit(2)

    # PROCESS THE @ARGS

    for opt, arg in opts:
        
        if opt == '-h': # help file                                                 \/-- @@@
        
            print 'collapseMotifTree_progressiveMode.py -s <species threshold> -H <entropy threshold>\n\nDescription:\n\nA wrapper for obtaining the blacklist at given species number threshold, S, and entropy threshold, H, (CTRL+F: blacklist_criteria_check) and e value (CTRL+F: @e-value)\n\nArguments:\n\n e = 0.05                      # CTRL+F: @e-value\n\nspecies_threshold = 3.0       # CTRL+F: blacklist_criteria_check\n\nentropy_threshold = 1.0       # CTRL+F: blacklist_criteria_check\n\nExample:\n\nblacklist_then_summaryStats_from_shell( e = 0.05, species_threshold=3.0, entropy_threshold=1.0)\n\n'

            sys.exit()

        elif opt in ("-f","--file1"):    # long and short args, see                  <-- @@@

            f = arg # don't forget to correct the type!!                             <-- @@@

        elif opt in ("-s", "--file2"):   #                                           <-- @@@

            s = arg #                                                                <-- @@@

    # CMD or PYTHON EXECUTE CONTROL: by default try run as if on cmd line, with input and output args ... else if it errors it may be running from within python, so then run default args: 

    # ... script was run from command line
    try:

        # SCRIPT
        diff(f,s)           # <-- main script, @REPLICATE exactly as below           <-- @@@
        #                            ^@args here processed from above

    # ... script was run within Python, e.g. %paste in IPython or imported
    except UnboundLocalError:
        
        f = raw_input("first file to compare? e.g. file1 = '/home/ab108/popbio/data_andy/isatabs/andy_0411_MAP/old/Europe/a_collection.txt' ")
        
        s = raw_input("second file to compare? e.g. file2 = '/home/ab108/0VB/isatab/malariaAtlas/data/isatab/Europe/a_collection.txt' ")

        diff(f,s)

if __name__ == "__main__":
    main(sys.argv[1:])

