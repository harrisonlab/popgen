#!/usr/bin/python
from sys import argv
from collections import defaultdict
#Takes in a dictionary of the form cultivar_id-> cultivar_name, and another in the form of sample_id -> cultivar_id, and renames samples in a new VCF file to append cultivar_name.

script, vcf, clone_list, sample_list = argv

out = vcf[:-4] + "_new_names.vcf"
out_h = open(out, 'w')

names = defaultdict()
with open (clone_list) as clone_list_h:
    for line in clone_list_h:
        lines = line.split()
        names[lines[0]] = lines[1]

samples = defaultdict()
with open (sample_list) as sample_list_h:
    for line in sample_list_h:
        lines = line.split()
        samples[lines[0]] = names[lines[1]]

all_names = []

with open (vcf) as vcf_h:
    for line in vcf_h:
        if line.startswith("##"):
            out_h.write(str(line))
        elif line.startswith("#"):
            fields = line.split()
            for f in fields[:9]:
                all_names.append(f)
            for f in fields [9:]:
                cultivar = samples[f]
                nn = f + "_" + cultivar
                all_names.append(nn)
            new_line = "\t".join(all_names)
            out_h.write(new_line + "\n")
        else:
            out_h.write(str(line))

out_h.close()