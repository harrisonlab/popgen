#! /usr/bin/env python

#Takes a file with a list of genes (1 gene/line) and a FASTA file.
#Outputs a FASTA file. Prints only sequences present in the list of genes
#First argument: list of genes
#Second argument: FASTA file
#Third arument: Reverse complement the sequence? Yes or No

from sys import argv
import os, sys, re
from Bio import SeqIO

script, gene_list, fasta_file, revcom = argv

genes = dict()

bare = r"(\w+)(.\w+$)"
filtered = r"\1.fasta"
o = re.sub(bare, filtered, gene_list)
out = open(o, 'w')

gene_list_h = open(gene_list)
for line in gene_list_h:
    gene_id = line.strip()
    genes[gene_id] = 1
gene_list_h.close()


for seq_record in SeqIO.parse(fasta_file, "fasta"):
    gene_id = str(seq_record.id)
    if gene_id in genes:
        if revcom == "No":
            out.write(">" + gene_id + "\n" + str(seq_record.seq) + "\n")
        elif revcom == "Yes":
            reverse = seq_record.seq.reverse_complement()
            out.write(">" + gene_id + "\n" + str(reverse) + "\n")
    else:
        pass

out.close()
