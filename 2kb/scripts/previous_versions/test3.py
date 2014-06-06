
import os 
home = os.path.expanduser('~')
os.chdir(home)

masking_char = 'n'

file_out_tmp = '0VB/2kb/data/meme_data/in/' + 'aedes_aegypti_tmp.bg2'
file_out     = '0VB/2kb/data/meme_data/in/' + 'aedes_aegypti.bg2'

# PRUNE OUT MASKING 'N's    <= meme prohibits 'N's in the -bfile, so we must prune em out
file_in_tmp = open(file_out_tmp,'r') # Open file for reading
file_out    = open(file_out,'w')     # Open file for writing

print '\tpruning illegal characters...'
while True:
    line = file_in_tmp.readline()

    print(line)

    if line == "":                    # break when finished
        break
    if not masking_char in line:

        print('detected n: no...  '+line)

        file_out.write(line)
    # else:
    #     print('detected n: yes... '+line)

file_in_tmp.close()
file_out.close()