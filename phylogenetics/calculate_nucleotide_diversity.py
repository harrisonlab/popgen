#! /usr/bin/env python
import os, sys
from sys import argv
import glob
import dendropy
from dendropy.calculate import popgenstat

script, regex = argv

##loop over all fasta files in a dir and calculate basic sequence diversity stats (Pi - nucleotide diversity, number of segregating sites,
##number of pairwise differences) and print to file
##OUTPUTS: excel_stats.txt (table with the results per gene), sequence_stats.txt (verbal summary of the results per gene)
f = open('sequence_stats.txt', 'w')
e = open('excel_stats.txt', 'w')
files = ["fileID"]
pis = ["Pi"]
segs = ["SegregatingSites"]
pdiffs = ["PairwiseDiffs"]
lists = [files, pis, segs, pdiffs]
for filename in glob.iglob(regex):
    seqs = dendropy.DnaCharacterMatrix.get(
    path= filename,
    schema="fasta")
    files.append(filename)
    pi = dendropy.calculate.popgenstat.nucleotide_diversity(seqs, ignore_uncertain=True)
    pis.append(pi)
    seg = dendropy.calculate.popgenstat.num_segregating_sites(seqs, ignore_uncertain=True)
    segs.append(seg)
    pdiff =dendropy.calculate.popgenstat.average_number_of_pairwise_differences(seqs, ignore_uncertain=True)
    pdiffs.append(pdiff)
    
    f.write("%s\n" % filename)
    f.write("Nucleotide diversity for sequences in the file is %s\n" % pi)
    f.write("Number of segregating sites for sequences in the file is %s\n" % seg)
    f.write("Average number of pairwise differences for sequences in the file is %s\n" % pdiff)

for l in lists:
    line = ""
    for a in l:
        line = line + str(a) + "\t"
    e.write(line + "\n")

f.close()
e.close()
