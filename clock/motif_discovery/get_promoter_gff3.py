#! /usr/bin/env python
import os, sys, re
from sys import argv

#Write out the CDS regions corresponding to the beginning
#of the coding region in Ensembl gff files
#And modify them so that they correspond to a given chunk of the
#putative promoter upstream of the start codon.

#Write out plus and minus strands to separate files.
#After FASTA sequence extraction need to reverse complement sequences
#in the file containing minus strand sequences.
/home/armita/prog/emboss/EMBOSS-4.0.0/bin/revseq

#Also, write a conversion table with output coordinates
#and protein name to cross-reference with the
#the results from GFF to FASTA script from bedtools
