#! /usr/bin/env python
from sys import argv
from collections import defaultdict as dd
import os, sys, re

#Split the OrthoMCL output file with the output orthgroups into the
#same format as OrthoFinder and count up the mumber of genes per orthogroup per species.

#Inputs:
#1st argument: .txt file output from OrthoMCL analysis
#2d argument: list of genome names used in OrthoMCL analysis and present in the
# .txt output. Seperate each genome name by double comma, e.g. AA134,,GB56

########################
# Sample usage
# python rearrange_orthomcl_results2.py All_Strains_orthogroups.txt  A4,,Bc1,,Bc16,,Bc23,,Nov27,,Nov5,,Nov71,,Nov77,,Nov9,,ONT3,,SCRP245_v2
#######################3

#Outputs:
#File with suffix "_rearranged.tsv" - contains orthomcl results file rearranged
#to contain genome in columns, orthogroups in rows and gene names in individual cells
#File with suffix "_counts.tsv" - the same as above but with gene count numbers instead
#of gene names in individual cells

script, orthomcl_file, genome = argv
#Obtain a list containing all the genome names used in the OrthoMCL analysis
genomes = genome.split(",,")
#Regex to specify output files
input_f = r"(\w+)(.txt$)"
output_f = r"\1_rearranged.tsv"
output_c = r"\1_counts.tsv"
sumfile = re.sub(input_f, output_f, orthomcl_file)
sumfile2 = re.sub(input_f, output_c, orthomcl_file)
sum_out = open(sumfile, 'w')
sum_out2 = open(sumfile2, 'w')

#Dictionary of lists
lr = dd(list)

#Regular expression used to obtain genome name id and gene ids from the input file
split_id_gene = r"(\w+)(\|)(.+)"

#Write file headers
sum_out.write("Orthogroup" + "\t" + "\t".join(genomes) + "\n")
sum_out2.write("Orthogroup" + "\t" + "\t".join(genomes) + "\n")

#open input file
orthomcl_h = open(orthomcl_file)
#Analyse the input file line-by-line
for line in orthomcl_h:
    fields = line.split()
    #Write the Orthogroup name to output
    sum_out.write(fields[0] + "\t")
    sum_out2.write(fields[0] + "\t")
    #Empty the lr dictionary. It will hold the current line's results
    lr.clear()
    #Iterate through fields[1] until last
    for record in fields[1:]:
        #Search for the regular expression to fish out the gene names for a given
        #genome in the current orthogroup
        result = re.search(split_id_gene, record)
        #Add the results to the temporary dictionary for the orthogroup (key - genome name)
        #value (gene names stored in a list)
        lr[result.group(1)].append(result.group(3))
    for g in genomes:
        #Write out the results for the current orthogroup (gene names seperated by comma)
        sum_out.write(",".join(lr[g]) + "\t")
        #Write out the results for the current orthogroup (gene counts)
        sum_out2.write(str(len(lr[g])) + "\t")
    sum_out.write("\n")
    sum_out2.write("\n")

sum_out.close()
sum_out2.close()
orthomcl_h.close()
