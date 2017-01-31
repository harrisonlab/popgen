#! /usr/bin/env python
import os, sys, re
from sys import argv

from Bio import SeqIO
#Takes in a FASTA file, and outputs a FASTA file without N containing sequences.

script, fasta_file = argv

bare = r"(\w+)(.fasta|.fa)"
out_sub = r"\1_noN.fasta"
out = re.sub(bare, out_sub, fasta_file)
out_h = open(out, 'w')

for seq_record in SeqIO.parse(fasta_file, "fasta"):
    seq = str(seq_record.seq).upper()
    seqid = str(seq_record.id)
    if 'N' in seq:
        pass
    else:
        out_h.write('>' + seqid + "\n")
        out_h.write(seq + "\n")

out_h.close()
