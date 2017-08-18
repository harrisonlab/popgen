#!/usr/bin/python

#To convert map-based haplotypes to the IMPUTE2 format needed for reference panel option,
#Use the orig_A_B file to get haplotypes and the shapeit input file (pmap file) to get other info. We also need a file with the list of individual ids (one ind per line) and a string with population name.

import sys, re
from sys import argv
import gzip
script, ab_file, pmap_file, list_of_inds, pop_name = argv

def convert(x):
    if x == "A":
        return "0"
    elif x == "B":
        return "1"

#Read in the A_B file - convert A to 0 and B to 1.
bare = r"(\w+)(.pmap)"
out_sub1 = r"\1.hap.gz"
out_hap = re.sub(bare, out_sub1, pmap_file)

fh_ab = open(ab_file)
hap_out = gzip.open(out_hap, 'w')

#Produce the Hap file.
snp_ids = list()
for line in fh_ab:
    tok = line.strip().split()
    snp_ids.append(tok[0])
    converted_01 = [convert(x) for x in tok[1:]]
    out_line = " ".join(converted_01)
    hap_out.write(out_line + "\n")
number_of_ids = len(converted_01) / 2
hap_out.close()
fh_ab.close()

#Read in the SNPs positions
snp_positions = dict()
with open(pmap_file) as fh_pmap:
    for line in fh_pmap:
        lines = line.strip().split()
        snp_positions[lines[1]] = lines[3]

#Produce the Legend file.
out_sub2 = r"\1.legend.gz"
out_legend = re.sub(bare, out_sub2, pmap_file)
legend_out = gzip.open(out_legend, 'w')
legend_out.write("id position a0 a1\n")
#Only write out SNPs which are in the target population to be phased file.
for snp in snp_ids:
    if snp in snp_positions:
        legend_out.write(snp + " " + snp_positions[snp] + " " + "A" + " " + "B" + "\n")
legend_out.close()

#Read in the file with the list of individual ids
inds_ids = list()
with open(list_of_inds) as fh_list_of_inds:
    for line in fh_list_of_inds:
        lines = line.strip().split()
        inds_ids.append(lines[0])

#Produce the Sample file
out_sub3 = r"\1.sample"
out_sample = re.sub(bare, out_sub3, pmap_file)
sample_out = open(out_sample, 'w')
#Create fake ids etc.ab_file
sample_out.write("sample population group sex\n")
for z in inds_ids:
    sample_out.write(z + " " + pop_name + " " + "ananassa" + " " + "1" + "\n")
sample_out.close() 