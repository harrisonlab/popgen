#! /usr/bin/env python
import os
import re
import sys
from collections import defaultdict as dd
from sys import argv
#Given a file with contig lengths and locations of the selected repeat (output from identify_telomere_repeats.py), quantify the number of repeats at xx bp off the ends of the contig.

script, contig_lengths, locations, cutoff_bp = argv

bare = r"(\w+)(.bed$)"
out_sub = r"\1"
o = re.sub(bare, out_sub, locations)
out_name = o + "_" + str(cutoff_bp) + ".txt"
out = open(out_name, 'w')

#Determine the end boundry for each contig:
end_boundary = dict()
cl_h = open(contig_lengths)
for line in cl_h:
    fields = line.split()
    contig_id = fields[0].strip()
    boundary = int(fields[1]) - int(cutoff_bp)
    end_boundary[contig_id] = boundary
cl_h.close
#Parse the results file and save results to a dictionary
counts_start = dd(int)
counts_end = dd(int)
fl_h = open(locations)
for line in fl_h:
    fields = line.split()
    contig_id = fields[0].strip()
    loc = int(fields[1])
    if loc < int(cutoff_bp):
        counts_start[contig_id] += 1
    elif loc > end_boundary[contig_id]:
        counts_end[contig_id] += 1
fl_h.close()

keylist = end_boundary.keys()
keylist.sort()
out.write("Used " + cutoff_bp + " bp to define contig ends" + "\n")
out.write("contig name" + "\t" + "number of repeats at 5' end" + "\t" + "number of repeats at 3' end" + "\n")
for key in keylist:
    out.write(key + "\t" + str(counts_start[key]) + "\t" + str(counts_end[key]) + "\n")
out.close()