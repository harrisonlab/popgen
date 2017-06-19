#! /usr/bin/env python

#Parse SnpEff annotated VCF file with only synonymous variants (produced by SnpSift)
#to retrieve only four-fold degenerate synonymous (silent) variants

from sys import argv
import os, sys, re

script, vcf = argv

bare = r"(\w+)(.vcf)"
one = r"\1_4fd.vcf"
o = re.sub(bare, one, vcf)
vcf_out = open(o, 'w')

aa_deg = ["Ala", "Arg", "Gly", "Leu", "Pro", "Ser", "Thr", "Val"]
vcf_h = open(vcf)
for line in vcf_h:
    if line.startswith("#"):
        vcf_out.write(str(line))
    else:
        fields = line.split()
        ann = fields[7].split("|")
        aa = r"p.([a-zA-Z]+)(\d+)([a-zA-Z]+)"
        match = re.search(aa, ann[10])
        if match:
            b = match.group(1)
            if b in aa_deg:
                pos = r"(\d+)/(\d+)"
                num = re.search(pos, ann[11])
                a = int(num.group(1))
                if a % 3 == 0:
                    vcf_out.write(str(line))
                else:
                    pass
            else:
                pass

vcf_out.close()
vcf_h.close()
