#!/usr/bin/python
import sys
import re
import gzip
from sys import argv
from collections import defaultdict as dd
#First argument: Rob's VCF file to be converted into ped file (only containing markes from a single). 
#Second argument is a pmap file which contains the list of markers we are interested in to be output.
#Third argument is a string with population name.

script, vcffile, pmapfile, pop_name = argv
vcf_data = dd(lambda: dd(list))
vcf_alleles = dd()

def read_in_pmap(p):
    snps = dict()
    #Obtain the list of markers to be printed out in the PED file.
    with open (p) as p_h:
        for line in p_h:
            lines = line.split()
            snps[lines[1]] = lines[3]
        return snps
def save_snp_diploid(snp_id, ref, alt, indiv, f):
    if (f == "." or f == ".|."):
        vcf_data[indiv][snp_id] = ['?', '?']
    else:
        g = f.split("|")
        vcf_data[snp_id][indiv] = g
        vcf_alleles[snp_id] = [ref, alt]
          
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

def print_output(vcf_data, vcf_alleles, snps):
    #print LEGEND
    legend_file.write("id position a0 a1\n")
    keyl = snps.keys()
    keyl.sort()
    for snp in keyl: 
       if snp in vcf_data:
        legend_file.write(snp + " " + snps[snp] + " " + vcf_alleles[snp][0] + " " + vcf_alleles[snp][1] + "\n")
    #print SAMPLE
    sample_file.write("sample population group sex\n")
    for i in inds:
        sample_file.write(i + " " + pop_name + " " + "ananassa" + " " + "1" + "\n")
    #print HAP
    for snp in keyl: 
       if snp in vcf_data:
           for i in inds:
               hap_file.write(" ".join(vcf_data[snp][i]) + " ")
           hap_file.write("\n")
snps = read_in_pmap(pmapfile)        
vcf_handling(vcffile)
legend_file = gzip.open(pmapfile[:-5] + ".legend.gz", 'w')
sample_file = open(pmapfile[:-5] + ".sample", 'w')
hap_file = gzip.open(pmapfile[:-5] + ".hap.gz", 'w')
print_output(vcf_data, vcf_alleles, snps)
legend_file.close()
sample_file.close()
hap_file.close()