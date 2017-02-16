#!/usr/bin/python

'''
Eliminate repeat entries featuring intervals not present in the assembly file
'''

from sys import argv
import os, sys, re
from Bio import SeqIO

script, assembly, bed_file = argv

bare = r"(\w+)(.txt)"
one = r"\1_filtered.bed"
o = re.sub(bare, one, bed_file)
out = open(o, 'w')

gene_lim = dict()

for seq_record in SeqIO.parse(assembly, "fasta"):
    gene_id = str(seq_record.id)
    gene_len = len(seq_record.seq)
    gene_lim[gene_id] = gene_len

with open(bed_file, 'r') as bed:
    bed_data = bed.readlines()
    for line in bed_data:
        fields = line.split()
        if fields[0] in gene_lim:
            #print ("Rep:", fields[0], fields[2])
            #print (fields[0], gene_lim[fields[0]])
            if (float(fields[2]) <= gene_lim[fields[0]] and float(fields[1]) <= gene_lim[fields[0]]):
                out.write(line)