#! /usr/bin/env python
from sys import argv
import os, sys, re

script, input_file = argv
bare = r"(\w+)(.txt)"
one = r"\1.rlk"
o = re.sub(bare, one, input_file)
out = open(o, 'w')

gene_list_h = open(input_file)
for line in gene_list_h:
    var = 0
    line = line.strip()
    fields = line.split("\t")
    if not re.search("^\.$", fields[1]):
        var += 1
    if not re.search("^\.$", fields[2]):
        var += 1
    if not re.search("^\.$", fields[4]):
        var += 1
    if not re.search("^\.$", fields[5]):
        var += 1
    if var > 2:
        out.write(line + "\n")
out.close()
gene_list_h.close()
