#! /usr/bin/env python
import glob
import os, sys, re
from collections import defaultdict

all_matched = defaultdict()

for filename in glob.iglob('*summary.txt'):
    with open(filename) as f:
        for line in f:
            fields = line.split()
            if fields[0] in all_matched:
                all_matched[fields[0]] += int(fields[1])
            else:
                all_matched[fields[0]] = int(fields[1])


out = open("all_motifs.txt", 'w')
for k in all_matched:
    out.write(k + "\t" + str(all_matched[k]) + "\n")

out.close()

