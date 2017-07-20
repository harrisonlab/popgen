#!/usr/bin/python

#
# find best sample to use for each clone
# usage: zcat genotypes.tsv.gz | find_best_sample.py | gzip > best_sample.tsv.gz
#

import sys
import gzip

inpdata = sys.stdin

counts = {}

header = inpdata.readline().strip().split('\t')
assert header == ['sample_id','file','clone_id','name','alias_id','marker_id','snp_id','pipeline_id','genotype']
for line in inpdata:
    tok = line.strip().split('\t')
    sample_id = int(tok[0])
    clone_id = int(tok[2])

    #count number of calls per sample, group by clone id
    if not clone_id in counts: counts[clone_id] = {}

    if not sample_id in counts[clone_id]:
        counts[clone_id][sample_id] = 1
    else:
        counts[clone_id][sample_id] += 1

#for each clone_id print the sample_id with (possibly equal) most calls
print 'clone_id\tsample_id\tncalls'
for clone_id in counts:
    sample_list = [ [sample_id,counts[clone_id][sample_id]] for sample_id in counts[clone_id] ]
    sample_list.sort(key=lambda x:x[1],reverse=True)
    print '%d\t%d\t%d'%(clone_id,sample_list[0][0],sample_list[0][1])
