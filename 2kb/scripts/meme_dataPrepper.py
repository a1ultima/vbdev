

# Need to ensure data is ready for MEME

file_in     = open('../data/meme_data/test/anopheles_gambiae_upstream.fasta')                   # Open file for reading
file_out    = open('../data/meme_data/test/anopheles_gambiae_upstream_memeready_all.fasta','w')

print 'reading file...'
while True:
    header  = file_in.readline()
    sequence= file_in.readline()
    if header == "":                    # break when finished
        break
    #line_split = line.split('\t')  # split line where '\t' occurs
    #print(sequence)
    if len(sequence) < 8:
        print('gotcha')
        continue
    else:
        file_out.write(header+sequence)

file_in.close()
file_out.close()