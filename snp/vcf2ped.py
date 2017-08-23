#!/usr/bin/python
import sys
import re
from sys import argv
from collections import defaultdict as dd
from random import shuffle

#First argument: Rob's VCF file to be converted into ped file (only containing markes from a single). 
#Second argument is a pmap file which contains the list of markers we are interested in to be output.
#Third argument is either of:
# random: randomize the allele order (i.e. mix the haplotypes)
# nochange: output the same haplotype order as in the VCF file.

script, vcffile, pmapfile, switch = argv
vcf_data = dd(lambda: dd(list))

def read_in_pmap(p):
    snps = list()
    #Obtain the list of markers to be printed out in the PED file.
    with open (p) as p_h:
        for line in p_h:
            lines = line.split()
            snps.append(lines[1])
        return snps

def convert_genotype(ref, alt, genotype):
    if genotype == "0":
        return ref
    elif genotype == "1":
        return alt

def save_snp_diploid(snp_id, ref, alt, indiv, f):
    if (f == "." or f == ".|."):
        vcf_data[indiv][snp_id] = [0, 0]
    else:
        g = f.split("|")
        if switch == "random":
            shuffle(g)
            geno = [convert_genotype(ref, alt, x) for x in g]
        else:
            geno = [convert_genotype(ref, alt, x) for x in g]
        vcf_data[indiv][snp_id] = geno
          
inds = []
def vcf_handling(vcf_file):
    vcf_h = open(vcf_file)
    for line in vcf_h:
        if line.startswith("##"):
            pass
        elif line.startswith("#"):
            fields = line.split()
            #Index all the individuals
            for f in fields[9:]:
                inds.append(f)
        else:
            fields = line.split()
            #Iterate over individual genotypes
            for idx, f in enumerate(fields[9:]):
            #Cross-reference to individuals
                indiv = inds[idx]
                save_snp_diploid(fields[2], fields[3], fields[4], indiv, f)


def print_ped(vcf_data, snps):
    #Go by individual
    for ind in vcf_data.keys():
        ped_file.write("fam" + " " + ind + " 0 0 1 0 ")
        #Iterate over genotypes
        for snp in snps:
            if snp in vcf_data[ind]:
                ped_file.write(" ".join(vcf_data[ind][snp]) + " ")
            else:
                ped_file.write("0 0 ")
        ped_file.write("\n")
snps = read_in_pmap(pmapfile)        
vcf_handling(vcffile)
ped_file = open(pmapfile[:-5] + ".ped", 'w')
print_ped(vcf_data, snps)
ped_file.close()

