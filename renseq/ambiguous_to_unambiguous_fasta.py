#!/usr/bin/python

'''
In a FASTA file, change two-fold ambiguous IUPAC nucleotide characters (other than N and three-fold) to unambiguous (ATGC) selected at random.
'''

from sys import argv
import os, sys, re
from Bio import SeqIO
import random

script, fasta = argv


ambig_chars = {'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'], 'W': ['A', 'T'], 'K': ['G', 'T'], 'M': ['A', 'C']}

def convert_ambiguous(seq):
    seq_l = list(seq)
    seq_out = []
    for nuc in seq_l:
        if (nuc == "A" or nuc == "C" or nuc == "T" or nuc == "G" or nuc == "N"):
            seq_out.append(nuc)
        else:
            #pick 0 or 1 at random
            allele = round(random.random())
            seq_out.append(ambig_chars[nuc][allele])
    return ''.join(seq_out)


bare = r"(\w+)(.\w+$)"
filtered = r"\1_unambiguous.fasta"
o = re.sub(bare, filtered, fasta)
out = open(o, 'w')

for seq_record in SeqIO.parse(fasta, "fasta"):
    gene_id = str(seq_record.id)
    gene_seq = str(seq_record.seq)
    gene_seq = convert_ambiguous(gene_seq.upper())
    out.write(">" + gene_id + "\n")
    out.write(gene_seq + "\n")