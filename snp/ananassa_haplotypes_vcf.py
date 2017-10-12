#!/usr/bin/python
import sys
import re
from sys import argv
from collections import defaultdict

#Transform the output of the table from ananassa_genotypes_db.py to VCF. Reference position of SNPs given in a GFF file.
script, table, gff = argv

out = table + ".vcf"
out_h = open(out, 'w')

table_dict = defaultdict(dict)
markers = defaultdict()
individuals = defaultdict()
start = 0 

with open (table) as table_h:
    headers = table_h.readline().strip().split("\t")
    fieldno = len(headers)
    for line in table_h:
        lines = line.strip().split("\t")
        markers[lines[5]] = lines[-1]
        individuals[lines[2]] = 1
        for x in (range(1,fieldno)):
            table_dict[lines[5]][lines[2]] = lines[-2]

#Open the GFF file and print the output VCF file
positions = defaultdict()
with open (gff) as gff_h:
    for line in gff_h:
        lines = line.split()
        if not re.match(r'^#', line):
            (prefix, snp_id) = lines[8].split('=')
            positions[snp_id] = [lines[0], lines[3], lines[4]]
list_inds = sorted(individuals.keys())
out_h.write("##fileformat=VCFv4.2" + "\n")
out_h.write("#CHROM" + "\t" + "POS" + "\t" + "ID" + "\t" + "REF" + "\t" + "ALT" + "\t" + "QUAL" + "\t" + "FILTER" + "\t" + "INFO" + "\t" + "FORMAT" + "\t" + "\t".join(list_inds) + "\n")
#Iterate over markers and individuals
for m in markers.keys():
    #Obtain the chromosome
    if m in positions:
        chrom = positions[m][0]
        pos = positions[m][1]
    else:
        chrom = "Unknown"
        pos = str(start + 1)
    alleles = markers[m].split("/")
    ref = alleles[0]
    if ref == "": ref = "."
    alt = alleles[1]
    if alt == "": alt = "."
    out_h.write(chrom + "\t" + pos + "\t" + m + "\t" + ref + "\t" + alt + "\t" + "." + "\t" + "." + "\t" + "." + "\t" + "GT:DP" )
    #Iterate over individuals
    for i in list_inds:
        if i in table_dict[m]:
            genotype = table_dict[m][i]
            if genotype == "B|B":
                gen = "1|1:2"
            elif genotype == "A|A":
                gen = "0|0:2"
            elif genotype == "A|B":
                gen = "0|1:2"
            elif genotype == "B|A":
                gen = "1|0:2"
            else:
                gen = "."
        else:
            gen = "."
        out_h.write("\t" + gen)
    out_h.write("\n")   
out_h.close()