#! /usr/bin/env python
from sys import argv
import os, sys, re

#Add a column with InterProScan annotations to the duplication analysis output.
script, annotation_table, duplication_analysis = argv

o = os.path.basename(duplication_analysis) + "_ann"
out = open(o, 'w')

ids = {}

ff_h = open(annotation_table)
for line in ff_h:
    fields = line.split("\t")
    gene_id = "Fus2_" + fields[0].strip()
    ids[gene_id] = fields[27].strip() + "\t" + fields[30].strip()
ff_h.close

sf_h = open(duplication_analysis)
for line in sf_h:
    fields = line.split("\t")
    if fields[0] in ids:
        out.write(line.strip() + "\t" + ids[fields[0]] + "\n")
    else:
        out.write(line.strip() + "\t" + "NA" + "\n")
sf_h.close
out.close
