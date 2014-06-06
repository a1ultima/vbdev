

import subprocess 

#subprocess.call(['dust moo.fasta 10 > cow.fasta']) 

#subprocess.call(['dust moo.fasta 10']) 

f_out = open('cat.fasta','w')

subprocess.call(['dust','moo.fasta','10'],stdout=f_out) 


