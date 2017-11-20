#! /usr/bin/env python
from sys import argv
import os, sys, re

script, model_fc_output, htseq_output = argv
sample_names = htseq_output.split("_")

bare = r"(\w+)(.out)"
one = r"\1_fc.out"
o = re.sub(bare, one, htseq_output)

fc_ht_out = open(o, 'w')

#Read in model file from FeatureCounts
fc_data = dict()
with open (model_fc_output) as fc_data_h:
    for line in fc_data_h:
        fields = line.strip().split("\t")
        fc_data[fields[0]] = fields[1:]

#Substitute reads counts from that obtained with HTSeq counting. 
seen = dict()
with open (htseq_output) as htseq_output_h:
    for line in htseq_output_h:
        fields = line.strip().split("\t")
        if fields[0] in fc_data:
            fc_data[fields[0]][5] = fields[1]
            seen[fields[0]] = 1

#Print the output to file.
seen_k = sorted(seen.keys())

fc_ht_out.write("Geneid" + "\t" + "Chr" + "\t" + "Start" + "\t"	+ "End" + "\t" + "Strand" + "\t" + "Length" +	"\t" + "_".join(sample_names[0:5]) + "\n")
for k in seen_k:
    fc_ht_out.write(k + "\t" + "\t".join(fc_data[k]) + "\n")

fc_ht_out.close()


