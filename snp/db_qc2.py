#!/usr/bin/python
import sys
import re
from sys import argv
from collections import defaultdict

script, alias = argv
#Print my list of multiform and secondaryform markers.
#multiform means "Affx name with more than one marker id associated", secondaryform mean "Affx name with more than one marker id associated, *this* form does not match any istraw35 marker"

multi35 = open("my_multiform35.txt", 'w')
sec35 = open("my_secondaryform35.txt", 'w')
multi90 = open("my_multiform90.txt", 'w')
sec90 = open("my_secondaryform90.txt", 'w')

istraw35_markers = defaultdict(list)
istraw35_marker_ids = dict()
istraw90_markers = defaultdict(list)
istraw90_marker_ids = dict()

with open (alias) as alias_h:
    headers = alias_h.readline().strip().split("\t")
    for line in alias_h:
        lines = line.strip().split("\t")
        affx = lines[3]
        if lines[1] == "istraw90":
            istraw90_marker_ids[lines[4]] = 1
            if affx in istraw90_markers:
                istraw90_markers[affx].append(lines[4])
            else:
                istraw90_markers[affx] = [lines[4]]
        elif lines[1] == "istraw35":
            istraw35_marker_ids[lines[4]] = 1
            if affx in istraw35_markers:
                istraw35_markers[affx].append(lines[4])
            else:
                istraw35_markers[affx] = [lines[4]]

for affx in istraw35_markers:
    if len(istraw35_markers[affx]) > 1:
        print len(istraw35_markers[affx])
        markers = istraw35_markers[affx]
        markers = list(set(markers))
        multi35.write("\n".join(markers) + "\n")
        for m in markers:
            if m not in istraw90_marker_ids:
                sec90.write(m + "\n")

for affx in istraw90_markers:
    if len(istraw90_markers[affx]) > 1:
        print len(istraw90_markers[affx])
        markers = istraw90_markers[affx]
        markers = list(set(markers))
        multi90.write("\n".join(markers) + "\n")
    #Check for secondary form.
        for m in markers:
            if m not in istraw35_marker_ids:
                sec35.write(m + "\n")

multi35.close()
sec35.close()
multi90.close()
sec90.close()