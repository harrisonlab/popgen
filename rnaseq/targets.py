#! /usr/bin/env python
from sys import argv
import os, sys, re, glob
from collections import defaultdict as dd

batches = dd(int)
input_files = glob.glob('*.out')
print "ID", "\t", "File", "\t", "Mite",	"\t", "Plant", "\t", "Time", "\t", "Batch", "\t", "Label"
my_id = 0
for my_file in input_files:
    categories = my_file.split("_")
    my_id += 1
    print my_id, "\t", 
    print my_file, "\t", 
    mite = categories[2]
    print mite, "\t",  
    stress = categories[1]
    print stress, "\t",  
    time = categories[3]
    print time.replace("h", ""), "\t", 
    conditions = mite + stress + time
    if conditions in batches:
        batches[conditions] += 1
        print batches[conditions], "\t", 
    else:
        batches[conditions] = 1
        print batches[conditions], "\t", 
    label = mite[0] + stress[0] + str(time.replace("h", ""))
    print label
