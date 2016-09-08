#! /usr/bin/env python

#Retain only one (longest) alternative transcript emanating from the gene
#Python2.7

from sys import argv
import os, sys, re

from Bio import SeqIO

script, fasta_file = argv

bare = r"(\w+)(.fasta)"
one = r"\1_one.fasta"
o = re.sub(bare, one, fasta_file)
out = open(o, 'w')

#Dictionary holding all sequences to be printed to file
d = dict()

for seq_record in SeqIO.parse(fasta_file, "fasta"):
    seq_upper = str(seq_record.seq).upper()
    #strip seq_record.id of the .t\d portion at the end
    seq_name = str(seq_record.id)
    (seq_name_gene, trail) = seq_name.split(".")
    #If another transcript from the gene is already in the dictionary
    if seq_name_gene in d:
    #Check if the old or new entry is longer
        print (seq_name_gene)
        length_old = len(d[seq_name_gene][1])
        length_new = len(seq_upper)
        #If so subsitute the new entry for old
        if length_new > length_old:
            d[seq_name_gene] = [seq_name, seq_upper]
        else:
            pass
    #If no transcript from the gene available in the dict, add the current one
    else:
        d[seq_name_gene] = [seq_name, seq_upper]

#Print the dictionary of sequences to FASTA file
for k in sorted(d.iterkeys()):
    out.write(">" + d[k][0] + "\n" + d[k][1] + "\n")

out.close()
