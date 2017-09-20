import sys
import re
from sys import argv
from collections import defaultdict as dd
from ped_functions import *

script, ped, phenotype_scores = argv

samples = dd()
with open (phenotype_scores) as samples_tb_h:
    for line in samples_tb_h:
        lines = line.strip().split()
        samples[lines[0]] = lines[1]

infile = open(ped)
for line in infile:
    lines = line.split()
    if lines[1] in samples:
        lines[5] = samples[lines[1]]
    else:
        pass
    for l in lines:
        print l,
    print "\n"