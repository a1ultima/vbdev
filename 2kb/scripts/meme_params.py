
import subprocess
from datetime import datetime

#command = 'meme anopheles_gambiae_upstream_memeready_25.fasta -dna -mod anr -revcomp -minw 10 -maxw 12 -nmotifs 25 -maxsize 50000000 -p 4 -o ../out/ag'.split(' ')
#meme = 'meme anopheles_gambiae_upstream_memeready_25.fasta -dna -mod oops -revcomp -w 10 -nmotifs 1 -maxsize 50000000 -p 4 -o ../out/ag'.split(' ')

import os 
home = os.path.expanduser('~')
os.chdir(home)

in_bfile = '0VB/2kb/data/meme_data/in/'     + 'anopheles_gambiae.bg2'
in_fasta = '0VB/2kb/data/meme_data/test/'   + 'anopheles_gambiae_upstream_memeready_100_dusted_c5.fasta'
out_meme = '0VB/2kb/data/meme_data/out/'    + 'ag_anr'

meme_anr = [    'meme', 
                in_fasta, 
                '-dna', 
                '-mod', 'anr', 
                '-revcomp', 
                '-w', '10', 
                '-nmotifs', '10', 
                '-maxsize', '50000000', 
                '-p', '4', 
                '-bfile', in_bfile,
                '-oc', out_meme
                ]

all_meme = [meme_anr]

for meme in all_meme:
    start   = subprocess.Popen(['date'], stdout=subprocess.PIPE).communicate()[0]
    subprocess.call(meme)
    end     = subprocess.Popen(['date'], stdout=subprocess.PIPE).communicate()[0]
    tdelta  = datetime.strptime(end.split(' ')[3], '%H:%M:%S') - datetime.strptime(start.split(' ')[3], '%H:%M:%S')
    print(tdelta)