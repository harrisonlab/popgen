#! /usr/bin/env python
import os, sys, re
from sys import argv
from collections import defaultdict as dd
#Script to compare haplotypes preserved in IMPUTE2 files to that of a VCF to look for any differences.

script, vcffile, hapfile = argv

def convert_to_nt(g, ref, alt):
    if g == "0":
        return ref
    elif g == "1":
        return alt

#Read in the marker names from legend
legend_data = dd(list)
legend_file = open(hapfile[:-3] + "legend")
legend_file.readline() #skip header
for i,line in enumerate(legend_file):
    lines = line.split()
    legend_data[i] = lines

#Read in the sample names from sample
sample_data = dd(str)
sample_file = open(hapfile[:-3] + "sample")
sample_file.readline() #skip header
for i,line in enumerate(sample_file):
    lines = line.split()
    my_id = lines[0].replace("_phased", "")
    sample_data[i] = my_id

#Read in the haplotypes from haplotype file.
merged_data = dd(lambda: dd(list))
hap_file = open(hapfile)
for i,line in enumerate(hap_file):
    lines = line.split()
    #number of haplotypes
    length_l = len(lines)
    #Iterate over haplotypes:
    ind_pos = 0
    for haplo in range(0,length_l,2):
        genotype = [lines[haplo], lines[haplo+1]]
        ref = legend_data[i][2]
        alt = legend_data[i][3]
        genotype_n = [convert_to_nt(g, ref, alt) for g in genotype]
        #Retrieve marker name
        mn = legend_data[i][0]
        #Retrie individual name
        nn = sample_data[ind_pos]
        merged_data[mn][nn] = genotype_n
        #Move onto the next individual
        ind_pos += 1
#print merged_data
#Read in the data from VCF file and compare:
inds = []

vcf_h = open(vcffile)
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
        marker = fields[2]
        ref = fields[3]
        alt = fields[4]
            #Iterate over individual genotypes
        for idx, f in enumerate(fields[9:]):
        #Cross-reference to individuals
            indiv = inds[idx]
            genotype = f.split("|")
            genotype_n = [convert_to_nt(g, ref, alt) for g in genotype]
            if marker in merged_data:
                if indiv in merged_data[marker]:
                    if merged_data[marker][indiv] != genotype_n:
                        print "Error: no match for haplotype for %s in %s " % (marker, indiv)
                        print merged_data[marker][indiv]
                        print genotype_n
