#!/usr/bin/python

import os, sys, re
from sys import argv
from collections import defaultdict as dd

#Input file 1 = Table with raw counts
#Input file 2 = Table with gene CDS lengths
#Output file = Table with FPKM values

script, counts, lengths = argv
o = counts + ".fpkm"
out_h = open(o, 'w')

gene_lengths = dd()
with open (lengths) as lengths_h:
    headers1 = lengths_h.readline().strip().split("\t")
    for line in lengths_h:
        lines = line.strip().split("\t")
        gene_lengths[lines[0]] = int(lines[1])

count_samples1 = dd(dict)
count_samples2 = dd(dict)
with open (counts) as counts_h:
    headers2 = counts_h.readline().strip().split("\t")
    for line in counts_h:
        lines = line.strip().split("\t")
        for index, value in enumerate(lines[1:]):
            name = headers2[index+1]
            count_samples1[name][lines[0]] = int(value)
            count_samples2[lines[0]][name] = int(value)
#Print FPKM table.
#Calculate mapped library size for each sample.
lib_size = dd()
for s in count_samples1:
    counter = 0
    for c in count_samples1[s]:
        counter += count_samples1[s][c]
    #Save to new dict
    lib_size[s] = counter

#Print header
out_h.write("Gene")
for s in count_samples1:
    out_h.write("\t" + s)
out_h.write("\n")

#Print main body of the table
for c in count_samples2:
    out_h.write(c)
    for s in count_samples2[c]:
        if c in gene_lengths:
            fpkm_value = float(count_samples2[c][s]) / float(gene_lengths[c]) * 1000 / float(lib_size[s]) * 1000000 
            out_h.write("\t" + str(fpkm_value))
        else:
            out_h.write("\t" + "NA")
    out_h.write("\n")

out_h.close