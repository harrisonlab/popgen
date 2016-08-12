#! /usr/bin/env python
import os, sys
from sys import argv
import glob
import dendropy
from dendropy.calculate import popgenstat

script, filename = argv

##loop over all fasta files in a dir and calculate basic sequence diversity stats and print to file
f = open('sequence_stats.txt', 'a')
e = open('excel_stats.txt', 'a')

seqs = dendropy.DnaCharacterMatrix.get(
    path= filename,
    schema="fasta")
files.append(filename)
pi = dendropy.calculate.popgenstat.nucleotide_diversity(seqs, ignore_uncertain=True)
seg = dendropy.calculate.popgenstat.num_segregating_sites(seqs, ignore_uncertain=True)
pdiff =dendropy.calculate.popgenstat.average_number_of_pairwise_differences(seqs, ignore_uncertain=True)

f.write("%s\n" % filename)
f.write("Nucleotide diversity for sequences in the file is %s\n" % pi)
f.write("Number of segregating sites for sequences in the file is %s\n" % seg)
f.write("Average number of pairwise differences for sequences in the file is %s\n" % pdiff)

lists = [filename, pis, seg, pdiff]

for l in lists:
    line = ""
    for a in l:
        line = line + str(a) + "\t"
    e.write(line + "\n")

f.close()
e.close()
