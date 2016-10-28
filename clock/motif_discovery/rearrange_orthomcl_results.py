#! /usr/bin/env python
from sys import argv
from collections import defaultdict as dd
import os, sys, re

#Split the OrthoMCL output file with the output orthgroups into the
#same format as OrthoFinder.

#Specify an array containing the genome names in the order desired.
genomes = ["BoC", "Fus2", "FuG", "FuS", "MaO", "NeD", "NsC", "PoA", "SoM", "TrR", "VeA", "VeD"]

script, orthomcl_file = argv

sum_out = open("orthomcl_output.csv", 'w')
sum_out.write("Orthogroup" + "\t" + "\t".join(genomes) + "\n")
split_id_gene = r"(\w+)(\|)(.+)"

orthomcl_h = open(orthomcl_file)
for line in orthomcl_h:
    fields = line.split()
    sum_out.write(fields[0] + "\t")
    lr = dd(list)
    lr.clear()
    for record in fields[1:]:
        result = re.search(split_id_gene, record)
        lr[result.group(1)].append(result.group(3))
    for g in genomes:
        sum_out.write(",".join(lr[g]))
        sum_out.write("\t")
    sum_out.write("\n")

sum_out.close()
orthomcl_h.close()
