#!/usr/bin/python

#
# filter out genotype call not from the best pipeline
# also drop null genotype calls
# usage: zcat genotypes.tsv.gz | filter_pipeline.py best_pipeline.gz | gzip > outfile.tsv.gz
#

import sys
import gzip

inpdata = sys.stdin
bestfile = sys.argv[1]

#load sample_id to best pipeline_id info
sid2best = {}
f = gzip.open(bestfile)
header = f.readline().strip().split('\t')
assert header == ['sample_id','best_id']
for line in f:
    tok = line.strip().split('\t')
    sample_id = int(tok[0])
    best_id = int(tok[1])
    sid2best[sample_id] = best_id
f.close()

#filter genotype data
header_line = inpdata.readline().strip()
header = header_line.split('\t')
assert header == ['sample_id','file','clone_id','name','alias_id','marker_id','snp_id','pipeline_id','genotype']
print header_line
for line in inpdata:
    tok = line.strip().split('\t')

    #filter out null genotype calls
    if tok[8] == 'NULL': continue

    #filter out non-best pipeline ids
    sample_id = int(tok[0])
    pipeline_id = int(tok[7])
    if sid2best[sample_id] != pipeline_id: continue

    print line.strip()
