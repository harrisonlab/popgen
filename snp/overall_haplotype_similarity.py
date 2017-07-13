#!/usr/bin/python
import sys
import re
from sys import argv
from collections import defaultdict as dd
from ped_functions import *

#count how many loci they have the same allele call, ignoring homozygous loci, then if the agreement is less than 50% this implies that they were in fact using the opposite convention (ie haplo pipeline0 == genotype pipeline1)

script, ped1, ped2, info_file1, info_file2 = argv

bare = r"(^[a-zA-Z0-9]+)(.*$)"
out_sub = r"\1_overall.stat"
out = re.sub(bare, out_sub, ped1)

#Read in the marker definitions file 1
info_file1_pos = read_in_marker_definition(info_file1)
#Read in the marker definition file 2
info_file2_pos = read_in_marker_definition(info_file2)

#Read in haplotypes file 1 
haplotype1 = read_in_haplotype(ped1, info_file1_pos)
#Read in haplotypes file 2 
haplotype2 = read_in_haplotype(ped2, info_file2_pos)

out_h = open(out, 'w')
#Write header
out_h.write("sample_name" + "\t" + "matching_haplo" + "\t" + "all_haplo" + "\t" + "perc_matching_haplo" + "\n")
#Count the matches per chromosome per sample
for sample_1 in haplotype1:
    number_matches = 0
    number_markers = 0
    #Check if the individual present in both files being compared
    if sample_1 in haplotype2:
        #Check if the marker present in both files being compared
        for marker_1 in haplotype1[sample_1]:
            if marker_1 in haplotype2[sample_1]:
                number_markers += 1
                #Check if match
                if haplotype1[sample_1][marker_1] == haplotype2[sample_1][marker_1]:
                    number_matches += 1
                else:
                    #print haplotype1[sample_1][marker_1]
                    #print haplotype2[sample_1][marker_1]
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