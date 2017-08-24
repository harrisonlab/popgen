#!/usr/bin/python

#
# convert shapeit phased VCF file into crosslink-style phased loc file
# to enable haplotyping accuracy to be determined. No missing data allowed. 
# Use the corresponding pmap file for a given chromosome to filter variants.
#

import sys

vcffile = sys.argv[1]  #vcf file split by chrom
chrom = sys.argv[2]  #chrom name

def convert(ref, alt, haplos):
    if haplos == "0":
        return ref
    elif haplos == "1":
        return alt

pmap = open(chrom + ".pmap")
snp_ids = list()
for line in pmap:
    lines = line.split()
    snp_ids.append(lines[1])

f = open(vcffile)
of = open(chrom + "_original_vcf", 'w')
for line in f:
    if line.startswith("#"):
        pass
    else:
        gen_array = list()
        fields = line.split()
        uid = fields[2]
        if uid in snp_ids:
            alt = fields[3]
            ref = fields[4]
            for idx in fields[9:]:
                haplo = idx.split("|")
                gen_array.extend(haplo)
            codes = [convert(ref, alt, x) for x in gen_array]
            of.write(uid + ' ' +  ' '.join(codes) + "\n")
f.close()
of.close()
