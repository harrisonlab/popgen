#!/usr/bin/python

#
# split markers into one file per LG, drop markers not in the map at all
# drop markers with low call rate
#

import sys
import gzip
from os.path import basename
from collections import defaultdict

mapfile = sys.argv[1]
snpfile = sys.argv[2]
outdir  = sys.argv[3]

#load snp data
f = gzip.open(snpfile)
header = f.readline().strip().split(',')
data = {}
for line in f:
    tok = line.strip().split(',')
    snp_id = tok[0].split('(')[0]
    data[snp_id] = tok[1:]
f.close()

uid2lg = {}

#used to convert lg name to an integer
conv = {'A':1,'B':2,'C':3,'D':4}

#split map into one file per lg in plink-style map format
prevlg = None
fout = None
fsnp = None
prevbp = -1000
adj=0
f = open(mapfile)
for line in f:
    uid,lg,bp = line.strip().split(',')

    if not uid in data: continue

    if lg != prevlg:
        if fout != None:
            fout.close()
            fsnp.close()

        fout = open(outdir+'/'+lg+'.pmap','wb')
        fsnp = open(outdir+'/'+lg+'.tsv','wb')
        fsnp.write('\t'.join(header)+'\n')
        prevlg = lg

    ilg = int(lg[0]) * 10 + conv[lg[1]] #eg convert 1A to 11, 2D into 24 etc
    bp = int(round(float(bp)))
    cm = bp / 3.2e5                     #create fake cm positions at 320 kb / cM

    #prevent exactly the same bp positions from being produced
    if bp == prevbp:
        adj += 1
        bp += adj
    else:
        prevbp = bp
        adj = 0

    fout.write("%d %s %.8f %d\n"%(ilg,uid,cm,bp))
    fsnp.write(uid + '\t' + '\t'.join(data[uid]) + '\n')

    uid2lg[uid] = [lg,int(bp)]
f.close()
fout.close()
