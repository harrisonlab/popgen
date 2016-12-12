#! /usr/bin/env python
from sys import argv
import os, sys, re
#Remove terminal codon from CDS sequences in the input FASTA file,
#if it matches a stop codon.
#Output: FASTA file with "_ns" suffix

from Bio import SeqIO
script, fasta_file = argv

bare = r"(\w+)(.fasta|.fa)"
one = r"\1_ns.fasta"
o = re.sub(bare, one, fasta_file)
out = open(o, 'w')

for seq_record in SeqIO.parse(fasta_file, "fasta"):
    seq_upper = str(seq_record.seq).upper()
    seq_name = str(seq_record.id)
    if seq_upper[-3:] == "TAG":
        out.write(">" + seq_name + "\n")
        out.write(seq_upper[:-3] + "\n")
    elif seq_upper[-3:] == "TAA":
        out.write(">" + seq_name + "\n")
        out.write(seq_upper[:-3] + "\n")
    elif seq_upper[-3:] == "TGA":
        out.write(">" + seq_name + "\n")
        out.write(seq_upper[:-3] + "\n")
    else:
        out.write(">" + seq_name + "\n")
        out.write(seq_upper + "\n")
out.close()
