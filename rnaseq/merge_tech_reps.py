#! /usr/bin/env python
from sys import argv
import os, sys, re
from collections import defaultdict as dd

script, sample_ass, my_suffix = argv

sample_dict = dd(list)
with open (sample_ass) as sample_ass_h:
    for line in sample_ass_h:
        fields = line.strip().split("\t")
        sample_dict[fields[0]].append(fields[1])

#print sample_dict
#Iterate over each biological sample:

for key in sample_dict:
    counts = dd(int)
    other_fields = {}
    for sample in sample_dict[key]:
        input_file = sample + my_suffix
        print input_file
        with open (input_file) as input_file_h: 
            header = input_file_h.readline().split("\t")
            for line in input_file_h:
                fields = line.strip().split("\t")
                counts[fields[0]] += int(fields[6])
                other_fields[fields[0]] = fields[1:6]
    old_sample_id = header[6].split("_")
    new_sample_id = old_sample_id[0:4]
    header[6] = "_".join(new_sample_id) + "_" + key
    write_out = open (header[6] + "_fc.out", 'w') 
    write_out.write("\t".join(header) + "\n")
    my_keys = sorted(counts.keys())
    for a in my_keys:
        write_out.write(a + "\t" + "\t".join(other_fields[a]) + "\t" + str(counts[a]) + "\n")
    write_out.close()

