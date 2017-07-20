#!/usr/bin/python

#
# filter out genotype calls not from the best sample for each clone
# usage: zcat genotypes.tsv.gz | filter_sample.py best_sample.gz | gzip > outfile.tsv.gz
#

import sys
import gzip

inpdata = sys.stdin
bestfile = sys.argv[1]

#load clone_id to best sample_id info
cid2best = {}
f = gzip.open(bestfile)
header = f.readline().strip().split('\t')
assert header == ['clone_id','sample_id','ncalls']
for line in f:
    tok = line.strip().split('\t')
    clone_id = int(tok[0])
    sample_id = int(tok[1])
    cid2best[clone_id] = sample_id
f.close()

#filter genotype data
header_line = inpdata.readline().strip()
header = header_line.split('\t')
assert header == ['sample_id','file','clone_id','name','alias_id','marker_id','snp_id','pipeline_id','genotype']
print header_line
for line in inpdata:
    tok = line.strip().split('\t')

    #filter out non-best sample ids
    sample_id = int(tok[0])
    clone_id = int(tok[2])
    if cid2best[clone_id] != sample_id: continue

    print line.strip()
