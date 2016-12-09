#! /usr/bin/env python
import os, sys, re
from sys import argv
#Process the output VCF file from fasta_to_vcf.sh script to convert the output VCF
#file to proper VCF format. Output file has "correct" suffix.

script, vcf_file = argv
vcf_h = open(vcf_file)

bare = r"(\w+)(.vcf)"
out_sub = r"\1_correct.vcf"
out = re.sub(bare, out_sub, vcf_file)
vcf_out = open(out, 'w')

gen=r"(\d/)(\d)(.*)"


for line in vcf_h:
    if line.startswith("#"):
        vcf_out.write(str(line))
    else:
        fields = line.split()
        for f in fields[:8]:
            vcf_out.write(f + "\t")
        vcf_out.write("GT:DP" + "\t")
        for f in fields[9:]:
            genotype = re.search(gen, f)
            vcf_out.write(genotype.group(2) + ":1\t")
        vcf_out.write("\n")


vcf_out.close()
vcf_h.close()
