#! /usr/bin/env python

#Takes a file with a list of genes (1 gene/line) and a FASTA file.
#Outputs a FASTA file. Prints only sequences NOT present in the list of genes
#First argument: list of genes
#Second argument: FASTA file

from sys import argv
import os, sys, re
from Bio import SeqIO

script, gene_list, fasta_file = argv

genes = dict()

bare = r"(\w+)(.fa$|.fasta$)"
filtered = r"\1_filtered.fasta"
o = re.sub(bare, filtered, fasta_file)
out = open(o, 'w')

gene_list_h = open(gene_list)
for line in gene_list_h:
    gene_id = line.strip()
    genes[gene_id] = 1
gene_list_h.close()


for seq_record in SeqIO.parse(fasta_file, "fasta"):
    gene_id = str(seq_record.id)
    flag = 0
    for f in genes:
        if re.search(f, gene_id):
            flag = 1
    if flag == 0:
        out.write(">" + gene_id + "\n" + str(seq_record.seq) + "\n")

out.close()
