#!/usr/bin/python
import sys
import re
from sys import argv
from collections import defaultdict
#Distribution of number of probes per the same Affx ID in istraw35 and istraw90.

script, alias = argv
dist35 = open("probes_distr_35.txt", 'w')
dist90 = open("probes_distr_90.txt", 'w')

istraw35_markers = defaultdict(list)
istraw90_markers = defaultdict(list)

with open (alias) as alias_h:
    headers = alias_h.readline().strip().split("\t")
    for line in alias_h:
        lines = line.strip().split("\t")
        affx = lines[3]
        if lines[1] == "istraw90":
            if affx in istraw90_markers:
                istraw90_markers[affx].append(lines[2])
            else:
                istraw90_markers[affx] = [lines[2]]
        elif lines[1] == "istraw35":
            if affx in istraw35_markers:
                istraw35_markers[affx].append(lines[2])
            else:
                istraw35_markers[affx] = [lines[2]]

#Print histogram values:
def print_output(my_dictionary, chip):
    temp_dict = defaultdict(int)
    for a in my_dictionary:
        number_of_probes = len(my_dictionary[a])
        temp_dict[number_of_probes] += 1
    if chip == "35":
        for v in temp_dict:
            dist35.write(str(v) + "\t" + str(temp_dict[v]) + "\n")
    elif chip == "90":
        for v in temp_dict:
            dist90.write(str(v) + "\t" + str(temp_dict[v]) + "\n")

print_output(istraw35_markers, '35')
print_output(istraw90_markers, '90')

dist35.close()
dist90.close()