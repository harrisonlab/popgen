#! /usr/bin/env python
import os, sys
from sys import argv
from collections import defaultdict
from Bio import SeqIO

##############################################################
#   Produce FASTA files with all individuals in the VCF file.
#   Each chromosome -> new FASTA file.
##############################################################


##############################################################
# Usage: python vcf_to_fasta.py <VCF FILE> <FASTA FILE>
# 1|2|3
# 1 - haploid input
# 2 - diploid input, output as two separate haplotypes for each ind.
# 3 - diploid input, output as one sequence with ambiguity codes for each ind.
##############################################################

#Input: VCF - haploid/diploid (phased or unphased, treated as phased, regardless)
#       Reference assembly FASTA file
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
            if (len(fields[3]) == len(fields[4])):
                #Iterate over individual genotypes
                for idx, f in enumerate(fields[9:]):
                    #Cross-reference to individuals
                    a = idx - 9
                    indiv = inds[a]
                    n = f.split(":")
                    #Pass to the subroutine: contig_name, position, ref allele, alt allele, size of the allele
                    #individual name, and their genotype
                    if switch == "1":
                        save_snp_haploid(fields[0], int(fields[1]), fields[3], fields[4], len(fields[3]), indiv, f)
                    elif (switch == "2" or switch == "3"):
                        save_snp_diploid(fields[0], int(fields[1]), fields[3], fields[4], len(fields[3]), indiv, f)
    vcf_h.close()

def save_snp_haploid(contig, pos, ref, alt, length, indiv, f):
    if f == ".":
        vcf_data[contig][pos][indiv] = 'N' * length
    else:
        n = f.split(":")
        #First check for 'hidden' missing genotypes
        if n[2] == ".":
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
        vcf_data[contig][pos][indiv] = [['N' * length], ['N' * length]]
    else:
        n = f.split(":")
        #First check for 'hidden' missing genotypes
        if n[2] == ".":
            Nstring = 'N' * length
            vcf_data[contig][pos][indiv] = Nstring + Nstring
        else:
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

def print_fasta_haploid():
    #Loop over chromosomes
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        #Dict to store contig sequence for each individual
        temp = defaultdict(str)
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
                    #Loop over all individuals and append the sequence
                    for my_i in inds:
                        temp[my_i] += vcf_data[seqid][n][my_i]
                    variant_length = len(vcf_data[seqid][n][my_i])
                    #Move down the reference sequence if complex variant encountered
                    skip_steps = variant_length - 1
            #Taking reference site, if no variant present
                else:
                    for my_i in inds:
                        temp[my_i] += my_seq[i]
                if skip_steps:
                    [iterable.next() for x in range(skip_steps)]
            # Print FASTA chromosome
            out = str(seq_record.id + ".fasta")
            out_h = open(out, 'w')
            for my_i in inds:
                out_h.write('>' + my_i + "\n")
                out_h.write(str(temp[my_i]) + "\n")
            out_h.close()
        else:
            #Skip monomorphic loci
            continue


def print_fasta_diploid2():
    #Loop over chromosomes
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        #Dict to store contig sequence for each individual
        temp = defaultdict(lambda: defaultdict(str))
        seq = str(seq_record.seq).upper()
        seqid = str(seq_record.id)
        #Iterate over the range of the contig sequence
        iterable = iter(range(len(seq)))
        #Turn the sequence into a dictionary for easy indexing
        my_seq = list(seq)
        #Iterate over the contig FASTA sequence
        if seqid in vcf_data:
            for i in iterable:
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
                            temp[my_i]['0'] += genotype[0:variant_length]
                            temp[my_i]['1'] += genotype[variant_length:]
                        else:
                            #First chromosome
                            temp[my_i]['0'] += genotype[0]
                            #Second chromosome
                            temp[my_i]['1'] += genotype[1]
                        #Move down the reference sequence if complex variant encountered
                else:
                    for my_i in inds:
                        temp[my_i]['0'] += my_seq[i]
                        temp[my_i]['1'] += my_seq[i]
                if skip_steps:
                    [iterable.next() for x in range(skip_steps)]
                            #Taking reference site, if no variant present
        # Print FASTA chromosome
            out = str(seq_record.id + ".fasta")
            out_h = open(out, 'w')
            for my_i in inds:
                #First Chromosome
                out_h.write('>' + my_i + "_1" + "\n")
                out_h.write(str(temp[my_i]['0']) + "\n")
                #Second chromosome
                out_h.write('>' + my_i + "_2" + "\n")
                out_h.write(str(temp[my_i]['1']) + "\n")
            out_h.close()
        else:
            #Skip monomorphic loci
            continue

if switch == "1":
    vcf_handling()
    print_fasta_haploid()
if switch == "2":
    vcf_handling()
    print_fasta_diploid2()
if switch == "3":
    vcf_handling()
    print_fasta_diploid1()
