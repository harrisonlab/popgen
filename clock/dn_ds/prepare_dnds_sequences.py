#! /usr/bin/env python
from sys import argv
from collections import defaultdict as dd
import os, sys, re

gene_ids = dd()
#Input: matrix with species_names and gene identifiers.
#Output: for each gene, two fasta files - one with just one (first) gene
#from each species; the second with all homologs from a given species.

script, matrix_file = argv
mf_h = open(matrix_file)
genes = mf_h.readline().strip().split()
for g in genes:
    no = genes.index(g)
    gene_ids[no] = g

for line in mf_h:
    fields = line.split("\t")
    for i in fields[0:]:
        flag = 0
        no2 = fields.index(i)
        output = gene_ids[no2] + '_all.fasta'
        output2 = gene_ids[no2] + '_one.fasta'
        temp = open(output, 'a')
        temp2 = open(output2, 'a')
        all_ids = [x.strip() for x in i.split(',')]
        for ident in all_ids:
            if ident:
                file_h = open(fields[0])
                it = iter(file_h)
                for line_x in it:
                    if re.search(ident.strip(), line_x):
                        a = line_x
                        b = it.next()
                        temp.write(a)
                        temp.write(b)
                        if flag == 0:
                            temp2.write(a)
                            temp2.write(b)
                        flag = 1
                file_h.close()
        temp.close()
        temp2.close()
mf_h.close()
