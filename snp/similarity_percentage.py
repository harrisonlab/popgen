#! /usr/bin/env python
import os, sys, re
from sys import argv
from collections import defaultdict

script, vcf_file = argv


vcf_h = open(vcf_file)

bare = r"(\w+)(.vcf)"
out_sub = r"\1_distance.log"
out = re.sub(bare, out_sub, vcf_file)

out_log = open(out, 'w')
names = []
matrix = []
snps = 0

for line in vcf_h:
    if line.startswith("##"):
        pass
    elif line.startswith("#"):
        fields = line.split("\t")
        for f in fields[9:]:
            names.append(f)
        w = len(names)
        matrix = [[0 for x in range(w)] for y in range(w)]

    else:
        snps += 1
        count_a = 0
        fields = line.split("\t")
        for f in fields[9:]:
            n = f.split(":")
            count_b = 0
            for h in fields[9:]:
                k = h.split(":")
                if k[0] == n[0]:
                    pass
                    matrix[count_a][count_b] += 1
                else:
                    pass
                count_b += 1
            count_a += 1
out_log.write('\t' + '\t'.join(names))
for r in range(w):
    out_log.write(names[r] + "\t")
    for c in range(w):
        distance = float(matrix[r][c]) / float(snps)
        out_log.write("%f" % distance + "\t")
        out_log.write("\n")
