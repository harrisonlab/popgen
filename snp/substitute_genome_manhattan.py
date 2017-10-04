#! /usr/bin/env python
from sys import argv

#Argument 1:
#Read in GFF file 
#Argument 2:
#A table (tab-delimited)
#Argument 3
#Column index where the SNP ID resides in the table
#Argument 4
#Column index where the chromosome resides in the table

##Will output a table with chromosome positions from the GFF file substituted for from the previous value. 

script, gff_file, my_table, col_snp, col_chrom = argv

col_snp = int(col_snp)
col_chrom = int(col_chrom)

gff_h = open(gff_file)

all_snps = dict()
gff_h = open(gff_file)
for line in gff_h:
    if not line.startswith("#"):
        fields = line.split("\t")
        chrom = str(fields[0])
        snp = fields[8].split("=")
        all_snps[snp[1].strip()] = chrom
gff_h.close()

#print all_snps
my_table_h = open(my_table)
for line in my_table_h:
    fields = line.split("\t")
    if fields[col_snp] in all_snps:
        fields[col_chrom] = all_snps[fields[col_snp]]
    print "\t".join(fields)

my_table_h.close()