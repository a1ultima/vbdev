import subprocess
from datetime import datetime
import os 

home = os.path.expanduser('~')
os.chdir(home)


program = 'dreme' 


species  = 'Anopheles gambiae'.lower().replace(' ','_')


in_fasta = '0VB/2kb/data/meme_data/test/'   + species + '_upstream_memeready_1000_dusted_c2.fasta'
out_meme = '0VB/2kb/data/meme_data/out/'    + species + '_anr'


if program == 'meme':
    in_bfile = '0VB/2kb/data/meme_data/in/'     + species + '.bg2'
    meme_anr = [    'meme',                     # MEME program
                    in_fasta,                   # input fasta data , e.g. 2kb upstream vb & dust masked
                    '-dna',                     # dna motifs
                    '-revcomp',                 # + and - strands
                    '-mod',     'anr',          # how many occurrences of that motif per sequence
                    '-w',       '10',           # motif size to scan for 
                    '-nmotifs', '10',           # no. motifs to output
                    '-maxsize', '50000000',     # max inputfile size in characters
                    '-evt',     '0.0001',        # if e-value of motif > -evt: skip this motif
                    '-p',       '5',            # no. cores
                    '-bfile',   in_bfile,       # background markov model, Bob's script generates these at order 2 
                    '-oc',      out_meme
                    ]
    all_meme = [meme_anr]

elif program=='dreme':

    #in_bfile = '0VB/2kb/data/meme_data/in/'     + species + '.bg2'

    meme_anr = [    'meme',                     # MEME program
                    in_fasta,                   # input fasta data , e.g. 2kb upstream vb & dust masked
                    '-dna',                     # dna motifs
                    '-revcomp',                 # + and - strands
                    '-mod',     'anr',          # how many occurrences of that motif per sequence
                    '-w',       '10',           # motif size to scan for 
                    '-nmotifs', '10',           # no. motifs to output
                    '-maxsize', '50000000',     # max inputfile size in characters
                    '-evt',     '0.0001',        # if e-value of motif > -evt: skip this motif
                    '-p',       '5',            # no. cores
                    '-bfile',   in_bfile,       # background markov model, Bob's script generates these at order 2 
                    '-oc',      out_meme
                    ]
    all_meme = [meme_anr]


for meme in all_meme:
    start   = subprocess.Popen(['date'], stdout=subprocess.PIPE).communicate()[0]
    subprocess.call(meme)
    end     = subprocess.Popen(['date'], stdout=subprocess.PIPE).communicate()[0]
    tdelta  = datetime.strptime(end.split(' ')[3], '%H:%M:%S') - datetime.strptime(start.split(' ')[3], '%H:%M:%S')
    print(tdelta)