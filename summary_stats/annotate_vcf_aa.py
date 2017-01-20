#! /usr/bin/env python
from sys import argv
import os, sys, re, argparse
#A script to annotate the AA (ancestral allele) field in the VCF file
#given genotypes of specified individuals contained in the file.

#Input:
#First argument: Input VCF file with variants
#Second argument: Specify ploidy of the organism (arguments accepted: 1 or 2)
#Third argument: A list of the outgroup sample names, seperated by double commas (example:
#A244,,B4521,,C456)
#Output:
#A VCF file with the AA field added. Suffix: "_vcf_aa.vcf"

#AA field is added only when all the alleles
#observed in the outgroup match that present in the focal group to reduce the number of
#misannotations. Missing data and sequence depth is not taken into account so
#the input VCF file has to be pre-filtered to retain only high confidence SNPs
#prior to the start of the analysis.

script, vcf, ploidy, samples = argv

outgroup=samples.split(",,")

bare = r"(\w+)(.vcf)"
out_sub = r"\1_vcf_aa.vcf"
out = re.sub(bare, out_sub, vcf)
out_h = open(out, 'w')
vcf_h = open(vcf)

outgroup_ind = []
counter_vcf = 0
counter_aa = 0

def write_output(fields, current_allele):
    for z in fields[:7]:
        out_h.write(z + "\t")
        #Check if previous annotation present in the field
    if fields[7] != ".":
        out_h.write(fields[7] + ";")      
    out_h.write("AA=" + ','.join(current_allele) + "\t")
    for z in fields[8:]:
        out_h.write(z + "\t")
    out_h.write("\n")

def allele_check(current_allele, g):
    for geno in g:
        if geno not in current_allele:
            current_allele.append(geno)

for line in vcf_h:
    if line.startswith("##"):
        out_h.write(line)
    elif line.startswith("#"):
        out_h.write(line)
        fields = line.split()
        #Index all the individuals
        for idx, f in enumerate(fields):
            if f in outgroup:
                outgroup_ind.append(idx)
            else:
                pass
    else:
        counter_vcf += 1
        fields = line.split()
        alleles = []
        current_allele = []
        g = []
        alleles = [fields[3], fields[4]]
        #Check that biallelic variant, otherwise discard:
        if any("," in al for al in alleles):
            pass
        else:
            #Find what the ancestral alleles are
            for i in outgroup_ind:
                #Check if no empty genotypes
                if re.match(".*\d.*", fields[i]):
                    n = fields[i].split(":")
                    #Extract genotypes
                    if ploidy == "1":
                        g.append(alleles[int(n[0])])
                        allele_check(current_allele, g)
                    elif ploidy == "2":
                        d = n[0].split("/")
                        g.extend([alleles[int(d[0])], alleles[int(d[1])]])
                        allele_check(current_allele, g)
            write_output(fields, current_allele)
            counter_aa += 1
vcf_h.close()
out_h.close()

print ("Out of %d variants in the file %s, %d were annotated with (an) ancestral allele(s)"
% (counter_vcf, vcf, counter_aa))
