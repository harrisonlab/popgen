#!/usr/bin/python
import sys
import re
from sys import argv
from collections import defaultdict as dd
from ped_functions import *

#chop the chromosome every 1% total length and calculate the percentage of haplotypes which agree 100%, then plot percentage agreement relative to the position in the genome.

script, ped1, ped2, info_file1, info_file2 = argv

bare = r"(^[a-zA-Z0-9]+)(.*$)"
out_sub = r"\1_matching_fragments.stat"
out = re.sub(bare, out_sub, ped1)

#Read in the marker definitions file 1
info_file1_pos = read_in_marker_definition_2(info_file1)
#Read in the marker definition file 2
info_file2_pos = read_in_marker_definition_2(info_file2)

#Calculate the interval length by which to chop the chromosome.
positions = info_file1_pos.values()
positions_sorted = sorted(positions, key=int)
interval = (int(positions_sorted[-1]) - int(positions_sorted[0])) / 20

#Read in haplotypes file 1
haplotype1 = read_in_haplotype(ped1, info_file1_pos)
#Read in haplotypes file 2
haplotype2 = read_in_haplotype(ped2, info_file2_pos)

def write_results(sample_1, interval, positions_sorted, number_matches, number_markers):
    #Calculate the number of matches per interval:
    k = int(positions_sorted[0]) + interval
    for i in range(k, int(positions_sorted[-1]), interval):
        start = i - interval
        end = i + 1
        current_range = range(start, end, 1)
        match_count = 0
        marker_count = 0
        for z in current_range:
            if str(z) in number_matches:
                match_count += 1
            if str(z) in number_markers:
                marker_count += 1
    #Print the results:
        if (marker_count == 0 and match_count == 0):
            percentage_ag = 'NA'
        else:
            percentage_ag = (float(match_count) / float(marker_count)) * 100
        out_h.write(sample_1 + "\t" + str(start) + "\t" + str(end) + "\t" + str(match_count) + "\t" + str(marker_count) + "\t" + str(percentage_ag) + "\n")

out_h = open(out, 'w')
#Write header
out_h.write("sample_name" + "\t" + "interval_start" + "\t" + "interval_end" + "\t" + "matching_haplo" + "\t" + "all_haplo" + "\t" + "perc_matching_haplo" + "\n")
#Count the matches
for sample_1 in haplotype1:
    number_matches = dd(int)
    number_markers = dd(int)
    if sample_1 in haplotype2:
        for marker_1 in haplotype1[sample_1]:
            if marker_1 in haplotype2[sample_1]:
                number_markers[marker_1] += 1
                #Check if match
                if haplotype1[sample_1][marker_1] == haplotype2[sample_1][marker_1]:
                    number_matches[marker_1] += 1
                else:
                    print ("Haplotypes of " + str(marker_1) + " do not agree")       
            else:
                print ("Marker " + str(marker_1) + " not found")
        #write results:
        write_results(sample_1, interval, positions_sorted, number_matches, number_markers)
    else:
        print ("Sample " + sample_1 + " not found")  

out_h.close()