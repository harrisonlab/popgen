#! /usr/bin/env python
import os, sys, re, argparse

ap = argparse.ArgumentParser(description='Filter VCF calls to remove "fake" SNPs with homozygous reference allele (0/0 or 0|0) genotype across all the genotypes in the file. Also, optionally remove SNPs with CHROM and POS loaded from an optional table.')
ap.add_argument('--i',required=True,type=str,help='VCF file to be filtered')
ap.add_argument('--o',required=True,type=str,help='Output filename for the filtered VCF file')
ap.add_argument('--e',required=False,type=str,help='File with two columns (contig and position) of the SNPs to be removed from the file')

args = ap.parse_args()

vcf_file=args.i
out=args.o
flag=0
if args.e:
    remove_table=args.e
    flag = 1
    to_remove = dict()
    with open (remove_table) as rt:
        for line in rt:
            lines = line.strip().split()
            target = "\t".join(lines[0:2])
            to_remove[target] = 1
def read_vcf():
    vcf_h = open(vcf_file)
    for line in vcf_h:
        if line.startswith("#"):
            vcf_out.write(str(line))
        else:
            keep = 0
            fields = line.split()
            for f in fields[9:]:
                #Cross-reference to individuals
                    f_a = f.split(":")
                    #Check there is a non-refernce allele somewhere:
                    if re.search('1', f_a[0]):
                        keep = 1
                    elif re.search('2', f_a[0]):
                        keep = 1
            if keep == 1:
                vcf_out.write(str(line))
            

def read_vcf_filter_sites():
    vcf_h = open(vcf_file)
    for line in vcf_h:
        if line.startswith("#"):
            vcf_out.write(str(line))
        else:
            keep = 0
            fields = line.split()
            #Check that the SNP is not in the table of the sites to be eliminated.
            target = "\t".join(fields[0:2])
            if target in to_remove:
                next
            else:
                for f in fields[9:]:
                #Cross-reference to individuals
                    f_a = f.split(":")
                    #Check there is a non-reference allele somewhere:
                    if digits.match(f_a[0]):
                        keep = 1
                if keep == 1:
                    vcf_out.write(str(line))

vcf_out = open(out, 'w')
if flag:
    read_vcf_filter_sites()
else:
    read_vcf()
vcf_out.close()