#! /usr/bin/env python
from sys import argv
import os, sys, re
from Bio import SeqIO

#Extracts protein sequences from a given orthogroup in OrthoFinder (or parsed orthomcl) output
#Saves one fasta per species (name: <gene_name>_<species_name>.fa).
#Input fasta files with whole genome protein sequences need to be in the current directory.
#INPUT:
#1st argument - file containing the OrthoFinder output for a Orthogroup in question (plus header)
#2nd argument - Orthogroup gene name (e.g. freq)

script, gene_list, gene_name = argv
gene_list_h = open(gene_list)
species = dict()
genes_names = dict()

header = gene_list_h.readline()
headers = header.split("\t")
for i, a in enumerate(headers):
    print (i, a)
    species[i] = a
    no_fields = len(species) - 1

for l in gene_list_h:
    line = l.strip()
    fields = line.split("\t")
    for index, b in enumerate(fields[1:no_fields], start=1):
        if len(b) > 0:
            genes = b.split(",")
            for x in genes:
                genes_names[x.strip()] = 1
            print(index)
            input_file = species[index]
            proteome = open(input_file)
            fasta_sequences = SeqIO.parse(open(input_file),'fasta')
            out = gene_name + "_" + input_file
            with open(out, "w") as f:
                for seq in fasta_sequences:
                    if seq.id in genes_names:
                        SeqIO.write([seq], f, "fasta")
                    else:
                        pass
