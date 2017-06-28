#!/usr/bin/python
import sys
import re
from sys import argv
from collections import defaultdict as dd

#take individual markers and record for how far along the chromosome they agree... plot agreement length relative to the SNP starting position in the genome.

script, ped1, ped2, info_file1, info_file2 = argv

bare = r"(^[a-zA-Z0-9]+)(.*$)"
out_sub = r"\1_match_length.stat"
out = re.sub(bare, out_sub, ped1)

def read_in_marker_definition (infile):
    count = 0
    dict_file = dd()
    with open (infile) as infile_h:
        for line in infile_h:
            lines = line.strip().split()
            dict_file[count] = lines[1]
            count += 1
    return dict_file

#Read in the marker definitions file 1
info_file1_pos = read_in_marker_definition(info_file1)
#Read in the marker definition file 2
info_file2_pos = read_in_marker_definition(info_file2)

positions = info_file1_pos.values()
positions_sorted = sorted(positions, key=int)

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

def write_results(sample_1, matches, nonmatches):
    all_matches = matches + nonmatches
    all_matches_sorted = sorted(all_matches, key=int)
    end = all_matches_sorted[-1]
    matches_sorted = sorted(matches, key=int)
    nonmatches_sorted = sorted(nonmatches, key=int)
    for m in matches_sorted:
        #Find non-matching markers further down the chromosome
        limit_matches = [i for i in nonmatches_sorted if i > m]
        if len(limit_matches) > 0:
            haplotype_end = limit_matches[0]
        else:
            haplotype_end = end
        haplo_len = float(haplotype_end - m) 
        max_possible_len = float (end - m) 
        #When reaching the end of chromosome
        if (max_possible_len == 0 and haplo_len == 0):
            max_possible_len = 1
            haplo_len = 1
        haplotype_len_relative = haplo_len / max_possible_len * 100
        out_h.write(sample_1 + "\t" + str(m) + "\t" + str(haplotype_end) + "\t" + str(haplotype_len_relative) + "\n")


out_h = open(out, 'w')
#Write header
out_h.write("sample_name" + "\t" + "marker_pos" + "\t" + "haplotype_end" + "\t" + "perc_max_possible_haplo_len" +  "\n")
#Count the matches
for sample_1 in haplotype1:
    matches = list()
    nonmatches = list()
    if sample_1 in haplotype2:
        for marker_1 in haplotype1[sample_1]:
            if marker_1 in haplotype2[sample_1]:
                pass
                #Check if match
                if haplotype1[sample_1][marker_1] == haplotype2[sample_1][marker_1]:
                    matches.append(int(marker_1))
                else:
                    nonmatches.append(int(marker_1))    
            else:
                print ("Marker " + str(marker_1) + " not found")
        #write results:
        write_results(sample_1, matches, nonmatches)
    else:
        print ("Sample " + sample_1 + " not found")  

out_h.close()