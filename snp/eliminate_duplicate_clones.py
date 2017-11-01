#!/usr/bin/python
import sys
import re
from sys import argv
from collections import defaultdict as dd

#The script takes in a file output by the script and outputs a thinned sample id list with only one sample id per clone (cultivar) - select sample with the highest number of genotypes. 

script, sample, table = argv

#Count the number of genotypes per sample and save count in a dictionary.
sample_count = dd(int)
with open (table) as table_h:
    headers = table_h.readline().strip().split()
    for line in table_h:
        lines = line.strip().split()
        sample_count[lines[2]] += 1

chosen_sample = dd(int)
nulls = list()
with open (sample) as sample_h:
    headers = sample_h.readline().strip().split()
    for line in sample_h:
        lines = line.strip().split()
        #For a given sample id, check if sample present in sample_count dictionary,
        #and if so, substitute if has a higher number of genotypes than any other available
        #sample from the same clone id.
        if lines[0] in sample_count:
            if lines[1] in chosen_sample:
                #Check for null clone id - if present automatically accept into analysis.
                if lines[1] == "NULL":
                    nulls.append(lines[0])
                elif sample_count[lines[0]] > sample_count[chosen_sample[lines[1]]]:
                    chosen_sample[lines[1]] = lines[0]
                else:
                    pass
            else:
                chosen_sample[lines[1]] = lines[0]

for a in chosen_sample.values():
    print a

for b in nulls:
    print b