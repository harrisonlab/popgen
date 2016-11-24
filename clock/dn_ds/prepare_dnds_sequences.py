#! /usr/bin/env python
from sys import argv
from collections import defaultdict as dd
import os, sys, re

sequences = dd()
#Input: matrix with species_names and gene identifiers.
#Output: for each gene, two fasta files - one with just one (first) gene
#from each species; the second with all homologs from a given species.

script, matrix_file = argv
mf_h = open(matrix_file)
genes = mf_h.readline().strip().split()
for g in genes:
    sequences[g] = 1

for line in mf_h:
    fields = line.split()
    print (fields[0])


mf_h.close()
