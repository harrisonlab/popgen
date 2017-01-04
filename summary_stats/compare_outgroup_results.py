#! /usr/bin/env python
from sys import argv
import os, sys, re, argparse
from collections import defaultdict as dd

#The script to parse AA fields from two independent annotations of a given genome
#using a) read mapping to the common reference genome (script annotate_vcf_aa.py) b) whole genome alignment
#of contigs from the focal species and sister species with mauve-parser (script annotate_gen_aa.py).

#Input:
#1st argument: first VCF file with AA field annotations (read mapping to common reference)
#2nd argument: second VCF file with AA field annotations  (whole genome alignment of contigs)
#3rd argument: the ploidy of the organism (arguments accepted: 1 or 2)
#4th argument: add the fake genotype entry (options: Y or N)


#Output:
#VCF file with consensus AA field annotation. Suffix: "_cons"

#The script will produce a final AA field annotation using the following heuristics:
#1) If AA fields congruent -> print the AA field.
#2) If two different AA identified -> remove the AA field
#3) If one AA identified using one method and two AA using the other method, including
#the allele in method one, then -> print tha AA field with two alleles.
#4) If AA field identified using only one method -> print the AA field

script, vcf1, vcf2, ploidy, fake = argv

bare = r"(\w+)(.vcf)"
out_sub = r"\1_cons.vcf"
out = re.sub(bare, out_sub, vcf1)
out_h = open(out, 'w')

mapping_results = dd(lambda: dd(list))
aa_match = r"(AA=)(.*)"
vcf_h = open(vcf1)

def empty_fake_genotype():
    if fake == "Y":
        out_h.write("." + "\t" + ".")

def write_no_aa(fields):
    for z in fields[:7]:
        out_h.write(z + "\t")
    out_h.write("." + "\t")
    for z in fields[8:]:
        out_h.write(z + "\t")
        empty_fake_genotype()
    out_h.write("\n")

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


#Read in AA in the first VCF file into memory
for line in vcf_h:
    if line.startswith("#"):
        pass
    else:
        fields = line.split()
        if re.match("AA=", fields[7]):
            al = re.search(aa_match, fields[7])
            alleles = sorted(al.group(2).split(","))
            mapping_results[fields[0]][fields[1]] = alleles

vcf_h2 = open(vcf2)
#Read in AA in the second VCF file and print results
for line in vcf_h2:
    if line.startswith("##"):
        out_h.write(line)
    elif line.startswith("#"):
        out_h.write(line.strip())
        if fake == "Y":
            out_h.write("\t" + "ancestral_1" + "\t" + "ancestral_2" + "\n")
    else:
        alleles_present = dict()
        fields = line.split()
        #Collect the alleles present:
        alleles_present[fields[3]] = "0"
        alleles_present[fields[4]] = "1"
        if re.match("AA=", fields[7]):
            al = re.search(aa_match, fields[7])
            alleles2 = sorted(al.group(2).split(","))
            alleles = mapping_results[fields[0]][fields[1]]
            #Check if AA defined in the previous VCF file
            if len(alleles) > 0:
                common = list(set(alleles).intersection(alleles2))
                #Check if more than 1 allele present in the difference
                diff = set(alleles).symmetric_difference(alleles2)
                #check if any common ancestral alleles shared. If yes, proceed.
                if len(common) > 0:
                    #Check if criterion 3 met for retaining AA met:
                    if len(diff) > 1:
                        #remove the AA field
                        write_no_aa(fields)
                    else:
                        #Union of two lists:
                        aa_allele = list(set(alleles + alleles2))
                        write_output(fields, aa_allele, alleles_present)
                else:
                    #remove the AA field
                    write_no_aa(fields)
            else:
                #Add the AA field only from VCF2
                out_h.write(line.strip())
                if fake == "Y":
                    add_genotype(ploidy, alleles2, alleles_present)
                out_h.write("\n")
        else:
            alleles = mapping_results[fields[0]][fields[1]]
            #Check if AA field present in VCF1 and if so add it.
            if len(alleles) > 0:
                write_output(fields, alleles, alleles_present)
            else:
                out_h.write(line.strip())
                empty_fake_genotype()
                out_h.write("\n")
