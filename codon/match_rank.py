#! /usr/bin/env python
from sys import argv
import os, sys, re

#Add column with gene expression rank to the output from codonw analysis

script, rank, codonw = argv
bare = r"(\w+)(.\w+$)"
e1 = r"\1"
e2 = e1
ff_o = re.sub(bare, e1, os.path.basename(rank))
sf_o = re.sub(bare, e2, os.path.basename(codonw))
o = ff_o + "_" + sf_o
out = open(o, 'w')

ids = {}

ff_h = open(rank)
for line in ff_h:
    fields = line.split()
    gene_id = "Fus2_" + fields[1].strip()
    ids[gene_id] = fields[0].strip()
ff_h.close

sf_h = open(codonw)
for line in sf_h:
    fields = line.split()
    match = re.split('\W+', fields[0], 1)
    if match[0] in ids:
        out.write(line.strip() + "\t" + ids[match[0]] + "\n")
    else:
        out.write(line.strip() + "\t" + "NA" + "\n")
sf_h.close
out.close
