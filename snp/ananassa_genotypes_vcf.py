#!/usr/bin/python
import sys
from sys import argv

#Transform the output of the table from ananassa_genotypes_db.py to VCF. Reference position of SNPs given in a GFF file.
script, table, gff = argv