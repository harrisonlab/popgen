#! /usr/bin/env python
from sys import argv
import os, sys, re, argparse

#The script to parse AA fields from two independent annotations of a given genome
#using a) read mapping to the common reference genome b) whole genome alignment
#of contigs from the focal species and sister species with mauve-parser.

#Input:
#1st argument: first VCF file with AA field annotations (read mapping to common reference)
#2nd argument: second VCF file with AA field annotations  (whole genome alignment of contigs)
#3rd argument: the ploidy of the organism (arguments accepted: 1 or 2)
#4th argument: add the fake genotype entry (options: Y or N)


#Output:
#VCF file with consensus AA field annotation. Suffix: "_cons"

#The script will produce a final AA field annotation using the following heuristics:
#If AA field congruent -> print the AA field.
#If two different AA identified -> remove the AA field
#If one AA identified using one method and two AA using the other method, including
#the allele in method one, then -> print tha AA field with two alleles.

#If AA field identified using only one method -> print the AA field

script, vcf1, vcf2, ploidy, fake = argv

bare = r"(\w+)(.vcf)"
out_sub = r"\1_cons.vcf"
out = re.sub(bare, out_sub, vcf1)
out_h = open(out, 'w')

mapping_results = dd(lambda: dd(list))
aa_match = r"(AA=)(.*)"
vcf_h = open(vcf1)
#Read in AA in the first VCF file into memory
for line in vcf_h:
    if line.startswith("#"):
        pass
    else:
        fields = line.split()
        if re.match("AA=", fields[7]):
            al = re.search(aa_match, fields[7])
            alleles = al.group(2).split(",")
            mapping_results[fields[0]][fields[1]] = alleles.sort()

#Read in AA in the second VCF file and print results
for line in vcf_h:
    if line.startswith("##"):
        out_h.write(line)
    elif line.startswith("#"):
        out_h.write(line.strip())
        if fake == "Y":
            out_h.write("\t" + "ancestral_1" + "\t" + "ancestral_2" + "\n")
    else:
        fields = line.split()
        if re.match("AA=", fields[7]):
            al = re.search(aa_match, fields[7])
            alleles2 = sorted(al.group(2).split(","))
            alleles = mapping_results[fields[0]][fields[1]]
            #Check if AA defined in the previous VCF file
            if len(alleles) > 0:
                common = list(set(alleles).intersection(alleles2))
                #check if any common ancestral alleles shared. If yes, proceed.
                if len(common) > 0:
                #Union of two lists:
                    c = set(alleles) + set(alleles2)
                else:
                    #remove the AA field
            else:
                #Add the AA field only from VCF2
        else:
            #Check if AA field present in VCF1 and if so add it.
