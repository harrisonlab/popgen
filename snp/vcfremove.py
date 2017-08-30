#!/usr/bin/python
import sys
from sys import argv

script, vcffile_A, vcffile_B = argv

#Takes VCF file A and B, and outputs only those sites with IDs found in file A from file B.
vcf_h = open(vcffile_A)
markers = set()
for line in vcf_h:
    if line.startswith("#"):
            pass
    else:
        fields = line.split()
        markers.add(fields[2])
vcf_h.close()


output = open(vcffile_B[:-4] + "_filtered.vcf", 'w')
with open(vcffile_B) as vcf_h2:
    for line in vcf_h2:
        if line.startswith("#"):
            output.write(line)
        else:
            fields = line.split()
            if fields[2] in markers:
                output.write(line)
output.close()

    