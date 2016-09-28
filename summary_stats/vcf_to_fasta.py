#! /usr/bin/env python
import os, sys, re, argparse
from Bio import SeqIO

##############################################################
#   Produce FASTA files with all individuals in the VCF file.
#   Each chromosome -> new FASTA file.
##############################################################

#Input: VCF - haploid/diploid (phased or unphased, treated as phased, regardless)
#       Reference assembly FASTA file
#Can handle complex SNPs. but only if ref and alt allele of the same length,
#so that it does not mess up the alignment. Can handle missing data.

script, vcf_file, fasta_file = argv
