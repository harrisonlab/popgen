#! /usr/bin/env python
import os, sys, re
from sys import argv
#In VCF files add cultivar name to VCF header containing sample ids. For reference, read in table with cultivar name in the first column and sample name in the second column.

script, table, vcf = argv

#Read in cultivar_name - sample_name mappings

my_mapping = {}

with open(table) as table_h:
    for line in table_h:
        lines = line.strip().split()
        my_mapping[lines[1]] = lines[0]

all_names = []
with open (vcf) as vcf_h:
    for line in vcf_h:
        if line.startswith("##"):
            print line.strip()
        elif line.startswith("#"):
            fields = line.strip().split()
            for f in fields[:9]:
                all_names.append(f)
            for f in fields [9:]:
                if f in my_mapping:
                    cultivar = my_mapping[f]
                    nn = str(f) + "_" + cultivar
                else:
                    nn = str(f)
                all_names.append(nn)
            new_line = "\t".join(all_names)
            print new_line
        else:
            print line.strip()
