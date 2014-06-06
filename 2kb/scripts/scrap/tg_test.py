
from tg_repeat_investigate import get_repeat_geneIds

moo = get_repeat_geneIds('../data/meme_data/in/anopheles_gambiae_upstream_memeready_all.fasta')

from cogent.db.ensembl import HostAccount, Genome, Species
account = HostAccount('localhost','vbuser','Savvas')
Species.amendSpecies('Anopheles gambiae','anopheles gambiae')
genome = Genome(Species='Anopheles gambiae',Release=73,account=account)
gene = genome.getGeneByStableId('AGAP004679')


locs =[genome.getGeneByStableId(g).Location for g in moo]


locdic = {}

for i in locs:
    chrome = i.CoordName
    strand = i.Strand
    start = i.Start
    end = i.End 
    if locdic.has_key(chrome):
        locdic[chrome].append((float(start)+float(end))/2)
    else:
        locdic[chrome] = [(float(start)+float(end))/2]


import matplotlib.pyplot as plt

hmm = locdic['X']

plt.hist(hmm)

plt.show()


