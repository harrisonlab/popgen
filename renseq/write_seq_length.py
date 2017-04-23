#!/usr/bin/python

#Process FASTA file and print two columns (comma delimited): sequence id, sequence length

import os, sys, re
from sys import argv
from Bio import SeqIO

script, fasta = argv

bare = r"(\w+)(.fa$|.fas$|.fasta$)"
one = r"\1_lengths.txt"
o = re.sub(bare, one, fasta)
out = open(o, 'w')

for seq_record in SeqIO.parse(fasta, "fasta"):
    gene_id = str(seq_record.id)
    gene_len = len(seq_record.seq)
    out.write(gene_id + "," + str(gene_len) + "\n")

out.close