#! /usr/bin/env python
import os, sys, re, argparse
from sys import argv
from collections import defaultdict as dd

script, vcf, snp_locations = argv

subgenomes = list()
marker_location = dd()

with open (snp_locations) as snp_locations_h:
    for line in snp_locations_h:
        lines = line.strip().split(",")
        a = lines[1] + ".vcf"
        if a not in subgenomes:
            subgenomes.append(a)
        marker_location[lines[0]] = lines[1]

filedata = {filename: open(filename, 'w') for filename in subgenomes}

#filedata['1A'].write("OK")

vcf_h = open(vcf)
for line in vcf_h:
    if line.startswith("#"):
        for file in filedata.values():
                file.write(str(line))
    else:
        fields = line.split("\t")
        if fields[2] in marker_location:
            gen_chrom = marker_location[fields[2]]
            gen_chrom_file = gen_chrom + ".vcf"
            filedata[gen_chrom_file].write((str(line)))
    
for file in filedata.values():
    file.close()
        