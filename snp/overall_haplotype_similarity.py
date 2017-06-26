#!/usr/bin/python
import sys
import re
from sys import argv
from collections import defaultdict as dd

#count how many loci they have the same allele call, ignoring homozygous loci, then if the agreement is less than 50% this implies that they were in fact using the opposite convention (ie haplo pipeline0 == genotype pipeline1)

script, ped1, ped2, info_file1, info_file2 = argv

bare = r"(^[a-zA-Z0-9]+)(.*$)"
out_sub = r"\1_overall.stat"
out = re.sub(bare, out_sub, ped1)

def read_in_marker_definition (infile):
    count = 0
    dict_file = dd()
    with open (infile) as infile_h:
        for line in infile_h:
            lines = line.strip().split()
            dict_file[count] = lines[0]
            count += 1
    return dict_file

#Read in the marker definitions file 1
info_file1_pos = read_in_marker_definition(info_file1)
#Read in the marker definition file 2
info_file2_pos = read_in_marker_definition(info_file2)


def read_in_haplotype (ped, info_file_pos):
    bare = r"(.*/)([a-zA-Z0-9-_]*.CEL$)"
    matched = r"\2"
    #Read in haplotypes file 1
    haplotype = dd(lambda: dd(str))
    with open (ped) as ped_h:
        for line in ped_h:
            count = 0
            lines = line.strip().split()
            sample_name =  re.sub(bare, matched, lines[1])
            for i,k in zip(lines[6::2], lines[7::2]):
                if i != k:
                    haplo = i + "|" + k
                    count += 1
                    haplotype[sample_name][info_file_pos[count]] = haplo
    return haplotype

#Read in haplotypes file 1
haplotype1 = read_in_haplotype(ped1, info_file1_pos)
#Read in haplotypes file 2
haplotype2 = read_in_haplotype(ped2, info_file2_pos)

out_h = open(out, 'w')
#Write header
out_h.write("sample_name" + "\t" + "matching_haplo" + "\t" + "all_haplo" + "\t" + "perc_matching_haplo" + "\n")
#Count the matches
for sample_1 in haplotype1:
    number_matches = 0
    number_markers = 0
    if sample_1 in haplotype2:
        for marker_1 in haplotype1[sample_1]:
            if marker_1 in haplotype2[sample_1]:
                number_markers += 1
                #Check if match
                if haplotype1[sample_1][marker_1] == haplotype2[sample_1][marker_1]:
                    number_matches += 1
                else:
                    print ("Haplotypes of " + str(marker_1) + " do not agree")       
            else:
                print ("Marker " + str(marker_1) + " not found")
        #print results:
     #Percentage matches:
        per_mat = (float(number_matches) / float(number_markers)) * 100
        out_h.write(sample_1 + "\t" + str(number_matches) + "\t" + str(number_markers) + "\t" + str(per_mat) + "\n")
    else:
        print ("Sample " + sample_1 + " not found")  

out_h.close()