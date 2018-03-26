#! /usr/bin/env python
import os, sys, re
from sys import argv

script, subset, table = argv

#First file: two-columned table containing subset of entries in the second file, and to be output from the second file. Match up by two first columns.
to_keep = dict()
with open (subset) as rt:
    for line in rt:
        lines = line.strip().split()
        target = "\t".join(lines[0:2])
        to_keep[target] = 1

with open (table) as t:
    for line in t:
        lines = line.strip().split()
        target = "\t".join(lines[0:2])
        if target in to_keep:
            print line.strip()