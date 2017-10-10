#!/usr/bin/python
import sys
import re
from sys import argv
from collections import defaultdict

#Take a given marker_id and print the number of genotypes available for a given marker_id on istraw35 and istraw90, and alias_id, Affx codes associated with it on each chip.

script, alias, genotype, sample = argv

alias_dict = defaultdict(list)
marker_dict = defaultdict(list)
with open (alias) as alias_h:
    headers = alias_h.readline().strip().split("\t")
    for line in alias_h:
        lines = line.strip().split("\t")
        #Index by alias id 
        alias_dict[lines[0]] = [lines[1], lines[3]]
        marker_dict[lines[4]].append(lines[0])
#Check whether sample genotyped with istraw35 or istraw90
sample_dict = defaultdict()
with open (sample) as sample_h:
    headers = sample_h.readline().strip().split("\t")
    for line in sample_h:
        lines = line.strip().split("\t")
        sample_dict[lines[0]] = lines[4]

genotype_dict =  defaultdict(lambda: defaultdict(int))
with open (genotype) as genotype_h:
    headers = genotype_h.readline().strip().split("\t")
    for line in genotype_h:
        lines = line.strip().split("\t")
        pipeline_id = lines[1]
        sample_id = lines[2]
        alias_id = lines[3]
        if pipeline_id == "1":
            genotype_dict[alias_id]["1"] += 1
        elif pipeline_id == "2":
            genotype_dict[alias_id]["2"] += 1

headers = ["marker_id", "istraw35_affx_codes", "istraw90_affx_codes", "istraw35_alias_ids",	"istraw90_alias_ids", "pipeline1_number_istraw35_genotypes", "pipeline2_number_istraw35_genotypes", "pipeline1_number_istraw90_genotypes", "pipeline2_number_istraw90_genotypes"]
print "\t".join(headers)
for markers in marker_dict:
    istraw35_aliases = list()
    istraw90_aliases = list()
    istraw35_codes = list()
    istraw90_codes = list()
    i35_counts_p1 = 0
    i90_counts_p1 = 0
    i35_counts_p2 = 0
    i90_counts_p2 = 0
    print markers, " \t",
    for aliases in marker_dict[markers]:
        if alias_dict[aliases][0] == "istraw90":
            istraw90_aliases.append(aliases)
            istraw90_codes.append(alias_dict[aliases][1])
        elif alias_dict[aliases][0] == "istraw35":
            istraw35_aliases.append(aliases)
            istraw35_codes.append(alias_dict[aliases][1])
        for i35 in istraw35_aliases:
            if i35 in genotype_dict:
                if "1" in genotype_dic[i35]:
                    i35_counts_p1 += genotype_dict[i35]["1"]
                if "2" in genotype_dict[i35]:
                    i35_counts_p2 += genotype_dict[i35]["2"]
        for i90 in istraw90_aliases:
            if i90 in genotype_dict:
                if "1" in genotype_dict[i90]:
                    i90_counts_p1 += genotype_dict[i90]["1"]
                if "2" in genotype_dict[i90]:
                    i90_counts_p2 += genotype_dict[i90]["2"]
    print ";".join(istraw35_codes), " \t" , ";".join(istraw90_codes), " \t", ";".join(istraw35_aliases), " \t", ";".join(istraw90_aliases), " \t",
    print i35_counts_p1, " \t", i35_counts_p2,  "\t", i90_counts_p1, " \t", i90_counts_p2