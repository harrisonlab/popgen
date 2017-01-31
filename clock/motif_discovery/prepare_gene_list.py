#! /usr/bin/env python

from sys import argv
import os, sys, re

script, gene_list, suffix = argv

gene_list_h = open(gene_list)

for line in gene_list_h:
    fields = line.split()
    out_name = str(fields[0]) + '.' + str(suffix)
    out = open(out_name, 'w')
    for f in fields[1:]:
        fs = f.split(',')
        out.write('\n'.join(fs))
        out.write('\n')
    out.close()
gene_list_h.close()
