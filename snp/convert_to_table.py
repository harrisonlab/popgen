#!/usr/bin/python

#
# convert database dump of genotypes into table format
# usage: convert_to_table.py genotype_file.gz > outtable.csv
#

import numpy as np
import sys
import gzip

genofile = sys.argv[1]

#collect distinct marker_ids and sample_ids
#and map sample_ids to file names and marker_ids to snp_ids
markers = {}
samples = {}

f = gzip.open(genofile)
header = f.readline().strip().split('\t')
assert header == ['sample_id','file','clone_id','name','alias_id','marker_id','snp_id','pipeline_id','genotype']
for line in f:
    tok = line.strip().split('\t')

    samples[int(tok[0])] = tok[1]
    markers[int(tok[5])] = tok[6]

f.close()

#sort and index marker and sample ids
marker_list = markers.keys()
marker_list.sort()
marker2pos = {x:i for i,x in enumerate(marker_list)}

sample_list = samples.keys()
sample_list.sort()
sample2pos = {x:i for i,x in enumerate(sample_list)}

#allocate an array for the genotype values, initialised to all -1 (missing)
nsamples = len(sample_list)
nmarkers = len(marker_list)
arr = -np.ones((nmarkers,nsamples),dtype=np.int8)

#how to convert from AA,AB,BB to numerical codes
conv =\
{
    'AA':0,
    'AB':1,
    'BA':1,
    'BB':2,
}

#read through file again, populating the array
f = gzip.open(genofile)
header = f.readline().strip().split('\t')
assert header == ['sample_id','file','clone_id','name','alias_id','marker_id','snp_id','pipeline_id','genotype']
for line in f:
    tok = line.strip().split('\t')
    row = marker2pos[int(tok[5])] #marker_id
    col = sample2pos[int(tok[0])] #sample_id
    genotype = conv[tok[8]]
    arr[row,col] = genotype
f.close()

#dump to stdout as csv
fout = sys.stdout
for x in sample_list: fout.write(',%s(%d)'%(samples[x],x))
fout.write('\n')
for i,x in enumerate(marker_list):
    fout.write('%s(%d),'%(markers[x],x))
    fout.write(','.join([str(y) for y in arr[i,:]]))
    fout.write('\n')
