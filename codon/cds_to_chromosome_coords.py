#! /usr/bin/env python

from sys import argv
import os, sys, re

from Bio import SeqIO

#Parse out genomic locations of all genes in a gff file and put in a table
#Output format
#GeneName   Chromosome  GeneStart   GeneEnd
script, gff_file, prefix = argv
# prefix = optional prefix to append to gene names when printing out the table

bare = r"(\w+)(.gff$|.gff3$)"
one = r"\1_gene_table.txt"
o = re.sub(bare, one, gff_file)
out = open(o, 'w')

gff_h = open(gff_file)

for line in gff_h:
    if not line.startswith("#"):
        fields = line.split("\t")
        chrom = str(fields[0])
        start = str(fields[3])
        stop = str(fields[4])
        if fields[2] == "gene":
            (ids, gene_name, rest) = re.split('\W+', fields[8], 2)
            out.write(prefix + gene_name.strip() + "\t" + chrom + "\t" + start + "\t" + stop + "\n")
        else:
            pass
gff_h.close
out.close
