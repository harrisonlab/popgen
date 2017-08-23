#!/usr/bin/python

#
# convert shapeit phased VCF file into crosslink-style phased loc file
# to enable haplotyping accuracy to be determined. No missing data allowed. Split by chromosome (assumes sorted file).
#

import sys

vcffile = sys.argv[1]  #shapedit phased hap file

def convert(haplos):
    if haplos == "0":
        return 'A'
    elif haplos == "1":
        return 'B'

f = open(vcffile)
open_fh = dict()
for line in f:
    if line.startswith("#"):
        pass
    else:
        gen_array = list()
        fields = line.split()
        chrom = fields[0]
        uid = fields[2]
        for idx in fields[9:]:
            haplo = idx.split("|")
            gen_array.extend(haplo)
        codes = [convert(x) for x in gen_array]
        if chrom not in open_fh:
            filename = chrom + "_original_vcf"
            open_fh[chrom] = open(filename, 'w')
        else:
            open_fh[chrom].write(uid + ' ' +  ' '.join(codes) + "\n")
f.close()
for ff in open_fh.values():
    ff.close()
