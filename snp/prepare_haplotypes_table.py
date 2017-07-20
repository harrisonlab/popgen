#!/usr/bin/python
import csv
from sys import argv
from collections import defaultdict as dd
import glob
import os, re

script, marker_info, sample_info, genotype_info, alias_info, markers_out = argv
markers = dd(list)
sample = dd(list)
alias = dd(list)
markers_elim = list()
#Load the marker names to be excluded.
with open (markers_out) as markers_out_h:
    for line in markers_out_h:
        lines = line.strip().split("\t")
        markers_elim.append(lines[0])
#Load the marker db info into a dictionary.
with open (marker_info) as marker_info_h:
    headers = marker_info_h.readline().strip().split("\t")
    for line in marker_info_h:
        lines = line.strip().split("\t")
        #Assign marker_id, probe_id, snp_id and alleles to snp id
        markers[lines[0]] = [lines[1], lines[2], lines[3]]
#Load the sample db info into a dictionary.
with open (sample_info) as sample_info_h:
    headers = sample_info_h.readline().strip().split("\t")
    for line in sample_info_h:
        lines = line.strip().split("\t")
        #Assign clone_id, file to sample id
        sample[lines[0]] = [lines[1], lines[2]]
with open (alias_info) as alias_info_h:
    headers = alias_info_h.readline().strip().split("\t")
    for line in alias_info_h:
        lines = line.strip().split("\t")
        #Assign probe_id, snp_id to alias id
        alias[lines[0]] = [lines[2], lines[3]]
#Load the individual haplotypes from all the chromosomes in the current folder.
bare = r"(.*/)([a-zA-Z0-9-_]*.CEL$)"
matched = r"\2"
out_h = open("haplotype_table.txt", 'w')
marker_names = dd()
haplotypes = dd(lambda: dd(str))
for filename in glob.glob('*.info'):
    markers_temp = list()
    chromosome = filename[:-5]
    filename_hap = os.path.splitext(filename)[0] + '.ped.phased'
    with open(filename, 'r') as f:
        for line in f:
            lines = line.strip().split()
            marker_names[lines[0]] = chromosome
            markers_temp.append(lines[0])
    #Open the corresponding file with haplotypes.
    with open(filename_hap, 'r') as f_hap:
        for line in f_hap:
            counter = 0
            lines = line.strip().split()
            sample_name =  re.sub(bare, matched, lines[1])
            for i,k in zip(lines[6::2], lines[7::2]):
                #Need to subsitute C with B 
                if i == "C":
                    i = i.replace("C", "B")
                if k == "C":
                    k = k.replace("C", "B")
                haplo = i + "|" + k
                marker_name = markers_temp[counter]
                counter = counter + 1
                haplotypes[sample_name][marker_name] = haplo
#Load genotype table and print output table
with open (genotype_info) as genotype_info_h:
    headers = genotype_info_h.readline().strip().split("\t")
    for line in genotype_info_h:
        lines = line.strip().split("\t")
        genotype_id, sample_id, alias_id, genotype = lines[0], lines[2], lines[3], lines[5]
        #Check if the marker is not in the list of markers to eliminated
        marker_name = alias[alias_id][1]
        if marker_name not in markers_elim:
        #Check if this genotype has been phased and alleles match.
            sample_name = sample[sample_id][1]
            if sample_name in haplotypes:
                if marker_name in haplotypes[sample_name]:
                    #Check that haplotypes match 
                    array1 = haplotypes[sample_name][marker_name].split("|")
                    array2 = list(genotype)
                    if sorted(array1) == sorted(array2):
                        out_h.write(genotype_id + "\t" + "3" + "\t")
                        chrom = marker_names[marker_name]
                        chroms = list(chrom)
                        out_h.write(chroms[0] + "\t" + chroms[1] + "\t")
                        out_h.write(haplotypes[sample_name][marker_name] + "\n")
        else:
            pass
        
out_h.close()