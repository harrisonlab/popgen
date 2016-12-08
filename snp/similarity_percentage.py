#! /usr/bin/env python
import os, sys, re
from sys import argv

# The only argument: VCF file to be used as input. Outputs a matrix with the percentage similarity of shared alleles between each sample,
# which is then to be visualised using the distance_matrix.R script.
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
        #Collect all the sample names
        for f in fields[9:]:
            names.append(f)
        #Get the dimensions of the matrix for frequency calculations
        w = len(names)
        #Create the matrix
        matrix = [[0 for x in range(w)] for y in range(w)]

    else:
        #Counter for the number of SNPs processed
        snps += 1
        #Counter to keep track of the individual in the x dimension
        count_a = 0
        fields = line.split("\t")
        for f in fields[9:]:
            n = f.split(":")
            #Counter to keep track of the individual in the y dimension
            count_b = 0
            for h in fields[9:]:
                k = h.split(":")
                #Compare SNP alleles and count the matches
                if k[0] == n[0]:
                    matrix[count_a][count_b] += 1
                else:
                    pass
                count_b += 1
            count_a += 1
#Write the matrix to file
out_log.write('\t' + '\t'.join(names))
for r in range(w):
    n = names[r].strip()
    out_log.write(n + "\t")
    for c in range(w):
        #Get the percentage similarity of shared alleles
        distance = float(matrix[r][c]) / float(snps) * 100
        out_log.write("%.3f" % distance + "\t")
    out_log.write("\n")
