#! /usr/bin/env python
import os, sys, re, argparse
from collections import defaultdict as dd
from Bio import SeqIO

#Prepare FASTA input for pairwise Ka/Ks comparisons. All pairwise comparisons for genes in a given orthogroups output. 
ap = argparse.ArgumentParser()
ap.add_argument('--o',required=True,type=str,help='OrthoFinder results table containing orthogroups of interest')
ap.add_argument('--col1',required=True,type=int,help='Column containing orthogroups for species 1')
ap.add_argument('--fasta1',required=True,type=str,help='FASTA file with CDS sequences for species 1')
ap.add_argument('--prefix1',required=True,type=str,help='Prefix to be added to sequences from species 1')
ap.add_argument('--col2',required=True,type=int,help='Column containing orthogroups for species 2')
ap.add_argument('--fasta2',required=True,type=str,help='FASTA file with CDS sequences for species 2')
ap.add_argument('--prefix2',required=True,type=str,help='Prefix to be added to sequences from species 2')

args = ap.parse_args()

#Open the file with orthogroup results
ortho_h = open(args.o)
ortho_groups = dd(lambda: dd(str))
all_ids_1 = list()
all_ids_2 = list()

for line in ortho_h:
    fields = line.split("\t")
    s1 = args.col1 - 1
    s2 = args.col2 - 1
    species1_id = fields[s1].split(",")
    #strip empty space
    species1_id = [x.strip() for x in species1_id]
    #remove empty elements
    species1_id = filter(None, species1_id)
    all_ids_1.extend(species1_id)
    species2_id = fields[s2].split(",")
    species2_id = [x.strip() for x in species2_id]
    species2_id = filter(None, species2_id)
    all_ids_2.extend(species2_id)
    #Check if both lists non-empty, ie. ortholog sequences available for both species
    if (len(species1_id) > 0 and len(species2_id) > 0):
        for id1 in species1_id:
            for id2 in species2_id:
                ortho_groups[fields[0]][id1] = id2


#Read in sequences for species 1 into memory 
sequences_species1 = dict()

for seq_record in SeqIO.parse(args.fasta1, "fasta"):
    gene_id = str(seq_record.id)
    if gene_id in all_ids_1:
        sequences_species1[gene_id] = str(seq_record.seq)

#Read in sequences for species 2 into memory 
sequences_species2 = dict()

for seq_record in SeqIO.parse(args.fasta2, "fasta"):
    gene_id = str(seq_record.id)
    if gene_id in all_ids_2:
        sequences_species2[gene_id] = str(seq_record.seq)

#Print out all possible pairwise comparisons to file.
ortho_names = ortho_groups.keys()
for n in ortho_names:
    keys_id = ortho_groups[n].keys()
    counter = 1
    for k in keys_id:
        open 