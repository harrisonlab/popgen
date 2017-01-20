#! /usr/bin/env python
import os, sys, re
from sys import argv
import random
from collections import defaultdict
from Bio import SeqIO

##############################################################
#   Produce FASTA files with all individuals in the VCF file.
#   Each chromosome -> new FASTA file.
##############################################################

##############################################################
# Usage: python vcf_to_fasta.py <VCF FILE> <FASTA FILE>
# <1|2|3>
# 1 - haploid input
# 2 - diploid input, output as two separate phased haplotypes for each ind.
# 3 - diploid input, output as one sequence with ambiguity codes for each ind.
##############################################################

#Input: 1) VCF - haploid/diploid (phased or unphased, treated as phased, regardless)
#       2) Reference assembly FASTA file
#Can handle complex SNPs. but only if ref and alt allele of the same length,
#so that it does not mess up the alignment. Can handle missing data.
#Output: FASTA file for each chromosome with the sequences for all individuals.
#For diploid sequences, two options: either two sequences per contig or one sequences
#with IUPAC ambiguity codes placed in.

script, vcf_file, fasta_file, switch = argv
#A list holding all the individual names in the VCF files
inds = []
#Dictionary of dictonaries of dictonaries to store SNP data
vcf_data = defaultdict(lambda: defaultdict(dict))

def multi(aa, aa2):
    alleles = [aa, aa2]
    #Check that biallelic variant, otherwise discard:
    if any("," in al for al in alleles):
        return 0
    else:
        return 1

def vcf_handling():
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
            #Check for equal length of alleles
            #Check if more than 2 alleles present, otherwise discard
            mgen = multi(fields[3], fields[4])
            #Check if variants of equal length and biallelic
            if (len(fields[3]) == len(fields[4]) and mgen == 1):
                #Iterate over individual genotypes
                for idx, f in enumerate(fields[9:]):
                    #Cross-reference to individuals
                    indiv = inds[idx]
                    n = f.split(":")
                    #Pass to the subroutine: contig_name, position, ref allele, alt allele, size of the allele
                    #individual name, and their genotype
                    if switch == "1":
                        save_snp_haploid(fields[0], int(fields[1]), fields[3], fields[4], len(fields[3]), indiv, f)
                    elif switch == "2":
                        save_snp_diploid(fields[0], int(fields[1]), fields[3], fields[4], len(fields[3]), indiv, f)
                    elif switch == "3":
                        save_snp_diploid_iupac(fields[0], int(fields[1]), fields[3], fields[4], len(fields[3]), indiv, f)
    print("Finished processing VCF file.")
    vcf_h.close()

def save_snp_haploid(contig, pos, ref, alt, length, indiv, f):
    if f == ".":
        vcf_data[contig][pos][indiv] = 'N' * length
    else:
        n = f.split(":")
        #First check for 'hidden' missing genotypes
        if (len(n) > 2 and n[2] == "."):
            vcf_data[contig][pos][indiv] = 'N' * length
        else:
            if n[0] == "0":
                vcf_data[contig][pos][indiv] = ref
            elif n[0] == "1":
                vcf_data[contig][pos][indiv] = alt
            else:
                vcf_data[contig][pos][indiv] = 'N' * length

def save_snp_diploid(contig, pos, ref, alt, length, indiv, f):
    if f == ".":
        Nstring = 'N' * length
        vcf_data[contig][pos][indiv] = Nstring + Nstring
    else:
        n = f.split(":")
        #First check for 'hidden' missing genotypes
        if (len(n) > 2 and n[2] == "."):
            Nstring = 'N' * length
            vcf_data[contig][pos][indiv] = Nstring + Nstring
        else:
            #if phased genotypes
            n[0] = re.sub('\|', '/', n[0])
            g = n[0].split("/")
            if (g[0] == "0" and g[1] == "0"):
                vcf_data[contig][pos][indiv] = ref + ref
            elif (g[0] == "1" and g[1] == "1"):
                vcf_data[contig][pos][indiv] = alt + alt
            elif (g[0] == "0" and g[1] == "1"):
                vcf_data[contig][pos][indiv] = ref + alt
            elif (g[0] == "1" and g[1] == "0"):
                vcf_data[contig][pos][indiv] = alt + ref
            else:
                vcf_data[contig][pos][indiv] = Nstring + Nstring

def save_snp_diploid_iupac(contig, pos, ref, alt, length, indiv, f):
    if f == ".":
        vcf_data[contig][pos][indiv] = 'N' * length
    else:
        n = f.split(":")
    #First check for 'hidden' missing genotypes
        if (len(n) > 2 and n[2] == "."):
            vcf_data[contig][pos][indiv] = 'N' * length
        else:
            #if phased genotypes
            n[0] = re.sub('\|', '/', n[0])
            g = n[0].split("/")
            if (g[0] == "0" and g[1] == "0"):
                vcf_data[contig][pos][indiv] = ref
            elif (g[0] == "1" and g[1] == "1"):
                vcf_data[contig][pos][indiv] = alt
            #This bit works only for simple mono-nucleotide SNPs
            #In other cases, simply select the allele at random
            elif ((g[0] == "0" and g[1] == "1") or (g[0] == "1" and g[1] == "0")):
                if ref == "A":
		            if alt == "C":
			            code = "M"
		            elif alt == "G":
			            code = "R"
		            elif alt == "T":
			            code = "W"
		            elif alt == "A":
			            code = "A"
                elif ref == "C":
		            if alt == "A":
			            code = "M"
		            elif alt == "G":
			            code = "S"
		            elif alt == "T":
			            code = "Y"
		            elif alt == "C":
			            code = "C"
                elif ref == "G":
		            if alt == "A":
			            code = "R"
		            elif alt == "C":
			            code = "S"
		            elif alt == "T":
			            code = "K"
		            elif alt == "G":
			            code = "G"
                elif ref == "T":
		            if alt == "A":
			            code = "W"
		            elif alt == "C":
			            code = "Y"
		            elif alt == "G":
			            code = "K"
		            elif alt == "T":
			            code = "T"
                else:
                    allele = round(random.random())
                    if allele == "0":
                        code = ref
                    elif allele == "1":
                        code = alt
                    else:
                        code = ref

                vcf_data[contig][pos][indiv] = code

            else:
                vcf_data[contig][pos][indiv] = 'N' * length

def print_fasta_haploid():
    #Loop over chromosomes
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        counter = 0
        #Dict to store contig sequence for each individual
        temp = defaultdict(list)
        seq = str(seq_record.seq).upper()
        seqid = str(seq_record.id)
        #Iterate over the range of the contig sequence
        iterable = iter(range(len(seq)))
        #Turn the sequence into a dictionary for easy indexing
        my_seq = list(seq)
        if seqid in vcf_data:
        #Iterate over the contig FASTA sequence
            for i in iterable:
            #convert to position in the sequence to cross-ref with dictionary
                n = i + 1
                skip_steps = 0
            #Check if site contains variant
                if n in vcf_data[seqid]:
                    #Loop over all individuals and append the sequence into dictionary
                    #Note: extending string terribly inefficient here.
                    #So first use dict to hold each sequence, and join at the end.
                    for my_i in inds:
                        temp[my_i].append(vcf_data[seqid][n][my_i])
                    variant_length = len(vcf_data[seqid][n][my_i])
                    #Move down the reference sequence if complex variant encountered
                    skip_steps = variant_length - 1
            #Taking reference site, if no variant present
                else:
                    for my_i in inds:
                        temp[my_i].append(my_seq[i])
                if skip_steps:
                    [iterable.next() for x in range(skip_steps)]
            # Print FASTA chromosome
            print ("Printing FASTA chromosome")
            out = str(seq_record.id + ".fasta")
            out_h = open(out, 'w')
            for my_i in inds:
                current_seq = ''.join(temp[my_i])
                out_h.write('>' + my_i + "\n")
                out_h.write(current_seq + "\n")
            out_h.close()
        else:
        #Skip monomorphic loci
            continue


def print_fasta_diploid():
    #Loop over chromosomes
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        counter = 0
        #Dict to store contig sequence for each individual
        temp = defaultdict(lambda: defaultdict(list))
        seq = str(seq_record.seq).upper()
        seqid = str(seq_record.id)
        #Iterate over the range of the contig sequence
        iterable = iter(range(len(seq)))
        #Turn the sequence into a dictionary for easy indexing
        my_seq = list(seq)
        #Iterate over the contig FASTA sequence
        if seqid in vcf_data:
            for i in iterable:
                counter += 1
                #print (counter)
                #convert to position in the sequence to cross-ref with dictionary
                n = i + 1
                skip_steps = 0
                #Check if site contains variant
                if n in vcf_data[seqid]:
                #Loop over all individuals and append the sequence
                        #First chromosome
                    for my_i in inds:
                        genotype = vcf_data[seqid][n][my_i]
                        variant_length = len(genotype) / 2
                        skip_steps = variant_length - 1
                        if variant_length > 1:
                            temp[my_i]['0'].append(genotype[0:variant_length])
                            temp[my_i]['1'].append(genotype[variant_length:])
                        else:
                            #First chromosome
                            temp[my_i]['0'].append(genotype[0])
                            #Second chromosome
                            temp[my_i]['1'].append(genotype[1])
                #Move down the reference sequence if complex variant encountered
                else:
                    for my_i in inds:
                        temp[my_i]['0'].append(my_seq[i])
                        temp[my_i]['1'].append(my_seq[i])
                if skip_steps:
                    [iterable.next() for x in range(skip_steps)]
                            #Taking reference site, if no variant present
        # Print FASTA chromosome
            print ("Printing FASTA chromosome")
            out = str(seq_record.id + ".fasta")
            out_h = open(out, 'w')
            for my_i in inds:
                #First Chromosome
                out_h.write('>' + my_i + "_1" + "\n")
                current_seq = ''.join(temp[my_i]['0'])
                out_h.write(current_seq + "\n")
                #Second chromosome
                out_h.write('>' + my_i + "_2" + "\n")
                current_seq = ''.join(temp[my_i]['1'])
                out_h.write(current_seq + "\n")
            out_h.close()
        else:
            #Skip monomorphic loci
            continue

if switch == "1":
    vcf_handling()
    print_fasta_haploid()
if switch == "2":
    vcf_handling()
    print_fasta_diploid()
if switch == "3":
    vcf_handling()
    print_fasta_haploid()
