#! /usr/bin/env python
import os, sys, re
from sys import argv
from collections import defaultdict as dd

#Substitute protein names into contigs obtained from extracting ranges with
#bedtools.

#Input arguments: 1 - fasta file, 2 - conversion table produced with get_promoter_gff3.py

script, fasta, coord = argv
fasta_h = open(fasta)
bare = r"(\w+)(.fasta)"
out = r"\1_n.fasta"
out_h = re.sub(bare, out, fasta)

out_fasta = open(out_h, 'w')
fasta_h = open(fasta)
coord_h = open(coord)

table = dd()

for line in coord_h:
    fields = line.split()
    table[fields[0]] = fields[1]

for line in fasta_h:
    if line.startswith(">"):
        name = line[1:-1]
        out_fasta.write(">" + table[name] + "\n")
    else:
        out_fasta.write(line)

out_fasta.close()
fasta_h.close()
coord_h.close()
