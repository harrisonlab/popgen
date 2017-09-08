#!/usr/bin/python
import sys
import re
from sys import argv
from collections import OrderedDict as odd
#Add sample id to the master strawberry table. Only perfect case-blind matches accepted. 
script, crown_rot_tb, samples_tb = argv

crown_rot = odd()
#Read in the list of the cultivars in the crown rot GWAS table.
with open (crown_rot_tb) as crown_rot_tb_h:
    for line in crown_rot_tb_h:
        fields = line.strip().split("\t")
        #Convert to lower-case:
        search_token = fields[0].lower()
        crown_rot[search_token] = line.strip()

#read in the samples table
samples = dd()
with open (samples_tb) as samples_tb_h:
    for line in samples_tb_h:
        lines = line.strip().split("\t")
        samples[line.lower()] = lines

#print crown_rot.keys()
#print the results
matching_ids = dd(list)
for my_id in crown_rot.keys():
    for my_sample in samples.keys():
        if re.search(my_id, my_sample):
            matching_ids[my_id].append(samples[my_sample][0])
        else:
            pass

for my_id in crown_rot.keys():
    if my_id in matching_ids:
        print crown_rot[my_id], "\t", ", ".join(matching_ids[my_id])
    else:
        print crown_rot[my_id]