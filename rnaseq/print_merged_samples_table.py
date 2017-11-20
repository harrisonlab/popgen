#! /usr/bin/env python
from sys import argv
import os, sys

from collections import defaultdict as dd

script, sample_ass = argv

sample_dict = dd(list)
with open (sample_ass) as sample_ass_h:
    for line in sample_ass_h:
        fields = line.strip().split("\t")
        sample_dict[fields[0]] = (fields[1])

print "Sample name", "\t", "Organism", "\t", "Stress", "\t",  "Mite", "\t", "Time", "\t", "Run"
#Iterate over each biological sample:
for key in sample_dict:
    sample_name = sample_dict[key].split("_")
    #print sample_name
    sample_id = "_".join(sample_name[0:4]) + "_" + key
    print sample_id, "\t", "\t".join(sample_name[0:4]), "\t", key


