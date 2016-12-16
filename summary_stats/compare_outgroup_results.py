#! /usr/bin/env python
from sys import argv
import os, sys, re, argparse

#The script to parse AA fields from two independent annotations of a given genome
#using a) read mapping to the common reference genome b) whole genome alignment
#of contigs from the focal species and sister species with mauve-parser.

#Input:
#1st argument: first VCF file with AA field annotations
#2nd argument: second VCF file with AA field annotations

#Output:
#VCF file with consensus AA field annotation. Suffix: "_cons" 

#The script will produce a final AA field annotation using the following heuristics:
#If AA field congruent -> print the AA field.
#If two different AA identified -> remove the AA field
#If one AA identified using one method and two AA using the other method, including
#the allele in method one, then -> print tha AA field with two alleles.
#If AA field identified using only one method -> print the AA field
