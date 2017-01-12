#! /usr/bin/env python
from sys import argv
import os, sys, re

script, input_file = argv
bare = r"(\w+)(.\w+)"
one = r"\1.strand"
o = re.sub(bare, one, input_file)
out = open(o, 'w')

gene_list_h = open(input_file)
for line in gene_list_h:
    if re.search("frame-", line):
        var = "-"
    elif re.search("frame+", line):
        var = "+"
    line = re.sub(r"_frame.*", "", line)
    out.write(line.strip() + "\t" + var + "\n")

out.close()
gene_list_h.close()
