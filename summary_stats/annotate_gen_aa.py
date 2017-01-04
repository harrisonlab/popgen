#! /usr/bin/env python
from sys import argv
import os, sys, re, argparse
from collections import defaultdict as dd

# The scripts adds the AA field and fake genotype entries (optional) obtained in the script mauve-parser.pl
# to the specified VCF file.

#Input:
#1st argument - output from mauve_parser.pl
#2nd argument - VCF input file with variants to be annotated
#3rd argument - ploidy of the organism
#4th argument - add the fake genotype entry (options: Y or N)

#Output:
#A VCF file with the AA field added. Suffix: "_gen_aa.vcf", as well as a fake genotype
#Fake ancestral genotypes will be added under "ancestral_1" and  "ancestral_2"

#Sample input
#AA will only be added if BOTH alleles present in the two outgroups (two
#last columns in the table with results) are present in the focal species.
#Otherwise, no ancestral allele will be added.

#A multiple genome alignment where the first sequence provided to Mauve acts as the reference.
 #YHet    291809  T       A       T
 #YHet    291810  T       T       T

script, mauve, vcf, ploidy, fake = argv

bare = r"(\w+)(.vcf)"
out_sub = r"\1_gen_aa.vcf"
out = re.sub(bare, out_sub, vcf)
out_h = open(out, 'w')
vcf_h = open(vcf)

counter_vcf = 0
counter_aa = 0

mauve_results = dd(lambda: dd(list))

def add_genotype(ploidy, current_allele, alleles):
    if ploidy == "1":
        if len(current_allele) == 1:
            out_h.write(alleles[current_allele[0]] + ":20" + "\t" + alleles[current_allele[0]] + ":20")
        elif len(current_allele) == 2:
            out_h.write(alleles[current_allele[0]] + ":20" + "\t" + alleles[current_allele[1]] + ":20")
    elif ploidy == "2":
        if len(current_allele) == 1:
            out_h.write(alleles[current_allele[0]] + "/" + alleles[current_allele[0]] + ":20")
            out_h.write("\t" + alleles[current_allele[0]] + "/" + alleles[current_allele[0]] + ":20")
        elif len(current_allele) == 2:
            out_h.write(alleles[current_allele[0]] + "/" + alleles[current_allele[1]] + ":20")
            out_h.write("\t" + alleles[current_allele[0]] + "/" + alleles[current_allele[1]] + ":20")

def write_output(fields, current_allele, alleles):
    for z in fields[:7]:
        out_h.write(z + "\t")
    out_h.write("AA=" + ','.join(current_allele) + "\t")
    for z in fields[8:]:
        out_h.write(z + "\t")
    if fake == "Y":
        add_genotype(ploidy, current_allele, alleles)
    out_h.write("\n")

#Read in mauve Output
mauve_h = open(mauve)
for line in mauve_h:
    fields = line.split()
    if (fields[3] in ['A', 'T', 'C', 'G'] and fields[4] in ['A', 'T', 'C', 'G']):
        mauve_results[fields[0]][fields[1]] = [fields[3], fields[4]]
mauve_h.close()

#Iterate over VCF file
for line in vcf_h:
    if line.startswith("##"):
        out_h.write(line)
    elif line.startswith("#"):
        out_h.write(line.strip())
        if fake == "Y":
            out_h.write("\t" + "ancestral_1" + "\t" + "ancestral_2")
        out_h.write("\n")
    else:
        switch = 0
        alleles_present = dict()
        counter_vcf += 1
        fields = line.split()
        alleles = [fields[3], fields[4]]
        #Check that biallelic variant, otherwise discard:
        if any("," in al for al in alleles):
            pass
        else:
            #Collect the alleles present:
            alleles_present[fields[3]] = "0"
            alleles_present[fields[4]] = "1"
            #Compare the alleles present here with that from mauve_h
            for a in mauve_results[fields[0]][fields[1]]:
                if a in alleles_present:
                    switch = 1
                else:
                    switch = 0
                    break
                    #If both the outgroup alleles present in the focal species, add AA field.
            if switch:
                if mauve_results[fields[0]][fields[1]][0] == mauve_results[fields[0]][fields[1]][1]:
                    aa_allele = [mauve_results[fields[0]][fields[1]][0]]
                else:
                    aa_allele = mauve_results[fields[0]][fields[1]]
                write_output(fields, aa_allele, alleles_present)
                counter_aa += 1
            else:
                if fake == "Y":
                    out_h.write(line.strip() + "\t" + "." + "\t" + "." + "\n")
                else:
                    out_h.write(line)
vcf_h.close()
out_h.close()

print ("Out of %d variants in the file %s, %d were annotated with (an) ancestral allele(s)"
% (counter_vcf, vcf, counter_aa))
