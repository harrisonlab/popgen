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

#Iterate over each biological sample:
for key in sample_dict:
    counts = dd(int)
    other_fields = {}
    for sample in sample_dict[key]:
        input_file = sample + my_suffix
        print input_file
        with open (input_file) as input_file_h: 
            for line in input_file_h:
                fields = line.strip().split("\t")
                counts[fields[0]] += int(fields[1])
    old_sample_id = input_file.split("_")
    new_sample_id = old_sample_id[0:4]
    header = "_".join(new_sample_id) + "_" + key
    write_out = open (header + ".out", 'w') 
    my_keys = sorted(counts.keys())
    for a in my_keys:
        write_out.write(a + "\t" + str(counts[a]) + "\n")
    write_out.close()

