#!/usr/bin/env python

import sys
import random

#weights=[.25,.25,.25,.25]
dna=['A','G','C','T']

def weighted_choice(weights,dna):
    totals = []
    running_total = 0

    for w in weights:
        running_total += w
        totals.append(running_total)

    rnd = random.random() * running_total
    for i, total in enumerate(totals):
        if rnd < total:
            return dna[i]

def dna_gen(length):
    seq=''
    for i in range(length):
        seq=seq+weighted_choice(weights,dna)
    return seq

def dna_gen2(reps,length,weights,dna):
    for i in range (reps):
        print (dna_gen(length))

if __name__=='__main__':
    reps    = int(sys.argv[1])
    length  = int(sys.argv[2])
    weights = [float (w) for w in sys.argv[3:6]]
    dna_gen2(reps,length,weights,dna)