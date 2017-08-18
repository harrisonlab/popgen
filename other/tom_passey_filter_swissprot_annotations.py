#! /usr/bin/env python
from sys import argv
from collections import defaultdict as dd
import os, sys, re

script, overlap_file, swissprot_file = argv
out_file = overlap_file + "_annot"

#Read in the annotation file
annot = dd()
f = open(swissprot_file)
for line in f:
    tok = line.strip().split("\t")
    annot[tok[0]] = [tok[13], tok[14]]
f.close()
#Read in the overlap file
of = open(overlap_file)
out_h = open(out_file, 'w')
for line in of:
    tok = line.strip().split("\t")
    ids = tok[-1]
    idds = ids.split(";")
    t_id = idds[0].replace("ID=", "")
    out_h.write(t_id + "\t" + annot[t_id][0] + "\t" + annot[t_id][1] + "\n") 
of.close()
out_h.close()
