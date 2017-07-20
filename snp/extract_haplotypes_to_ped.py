#!/usr/bin/python

#
# extract haplotypes to ped format
#

import sys
import numpy as np

pedfile = sys.argv[1]  #ped file
hapfile = sys.argv[2]  #haplotype file from shapeit
outfile = sys.argv[3]  #output file

lgname = pedfile.split('/')[-1].replace('.ped','')
#if lgname == '3A.2': exit()

#load all haplotype information
#one marker per line
conv = {'1':'A','2':'B','3':'C','4':'D'}
allele = ['-','-']
data = []
f = open(hapfile)
for line in f:
    tok = line.strip().split()
    lg,snpid,pos = tok[0:3]
    pos = int(pos)

    allele[0],allele[1] = tok[3],tok[4]

    if allele[0] == 'B': allele[0] = 'C'
    if allele[1] == 'B': allele[1] = 'C'

    calls = [allele[int(x)] for x in tok[5:]]
    #assert len(calls) == 2 * len(sampleids)

    data.append([lg[0]+conv[lg[1]],snpid,pos,calls])

f.close()

#process all ped file lines
#one sample per line
f = open(pedfile)
fout = open(outfile,'wb')
for i,line in enumerate(f):
    tok = line.strip().split()
    fam = tok[0]
    sampid = tok[1]
    aux = tok[2:6]
    pedcalls = tok[6:]

    out = []
    for j,row in enumerate(data):
        if pedcalls[2*j] == '0':
            #data missing
            out.append('0')
            out.append('0')
        else:
            #phased haplotype
            calls = row[3]
            out.append(calls[2*i])
            out.append(calls[2*i+1])

    fout.write(fam + ' ' + sampid + ' ' + ' '.join(aux) + ' ' + ' '.join(out) + '\n')
fout.close()
f.close()
