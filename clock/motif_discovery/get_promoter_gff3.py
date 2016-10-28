#! /usr/bin/env python
import os, sys, re
from sys import argv

#Write out the CDS regions corresponding to the beginning
#of the coding region in Ensembl gff files
#And modify them so that they correspond to a given chunk of the
#putative promoter upstream of the start codon.

#Also, write a conversion table with output coordinates
#and protein name to cross-reference with the
#the results from GFF to FASTA script from bedtools
