#! /usr/bin/env python
from sys import argv
import os, sys, re
from collections import defaultdict as dd

script, rbb, annotation = argv

rbb_matches = dict()
with open (rbb) as rbb_h:
    for line in rbb_h:
        fields = line.strip().split()
        #print fields[1]
        rbb_matches[fields[1]] = fields[0]
        
#print rbb_matches

with open (annotation) as annotation_h:
    for line in annotation_h:
        fields = line.strip().split("\t")
        if fields[0].strip() in rbb_matches:
            print rbb_matches[fields[0].strip()], "\t", fields[1]
        else:
            pass