#! /usr/bin/env python
import os, sys
from sys import argv

#First argument: BLAST output for DAGChain conversion
#Second argument: output of gff genome annotation parsing from cds_to_chromosome_coords.py

#Prepare input for DAGChainer using BLAST output from run_blast_reciprocal.pl

#DAGChainer input format
#query chromosome <tab> query gene name <tab> query hit 5' end <tab> query hit 3' end <tab>
#subject chromosome <tab> subject gene name <tab> subject hit 5' end <tab> subject hit 3' end <tab>
#P-value\n

script, input_file, chrom_coords = argv
output = input_file + "_dagchainer"
gff_table = dict()

chrom_h = open(chrom_coords)
#Save the GFF table results to dictionary
for line in chrom_h:
    fields = line.split("\t")
    gff_table[fields[0]] = [fields[1], fields[2], fields[3]]
chrom_h.close()

out = open(output, 'w')

#Skip first line (with annotation)
input_h = open(input_file)
next(input_h)
for line in input_h:
    fields = line.split("\t")
    #look up query chromosome
    (q_name_gene, trail) = fields[0].split(".")
    chrom = gff_table[q_name_gene][0]
    out.write(chrom + "\t")
    #query gene name
    out.write(fields[0] + "\t")
    #query hit 5' end
    qs = int(gff_table[q_name_gene][2]) + int(fields[6])
    out.write(str(qs) + "\t")
    #query hit 3' end
    qe = int(gff_table[q_name_gene][2]) + int(fields[7])
    out.write(str(qe) + "\t")
    #look up subject chromosome
    (s_name_gene, trail) = fields[1].split(".")
    chrom = gff_table[s_name_gene][0]
    out.write(chrom + "\t")
    #subject gene name
    out.write(fields[1] + "\t")
    strand = fields[17].strip()
    if strand == "plus":
    #subject hit 5' end
        ss = int(gff_table[s_name_gene][2]) + int(fields[8])
        out.write(str(ss) + "\t")
    #subject hit 3' end
        se = int(gff_table[s_name_gene][2]) + int(fields[9])
        out.write(str(se) + "\t")
    elif strand == "minus":
    #subject hit 5' end
        ss = int(gff_table[s_name_gene][2]) + int(fields[9])
        out.write(str(ss) + "\t")
    #subject hit 3' end
        se = int(gff_table[s_name_gene][2]) + int(fields[8])
        out.write(str(se) + "\t")
    #else:
    #    quit()
    #P-value
    pvalue = str(fields[10])
    out.write(pvalue + "\n")

out.close()
input_h.close()
