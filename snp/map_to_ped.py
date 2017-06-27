#!/usr/bin/python
import sys
import re
from sys import argv
from collections import defaultdict as dd

#First argument: map based haplotypes 
#Second argument: sample_id and filename obtained from the sample table.abs
#Third argument: mapping of marker names to positions on the chromosome. 

script, map_file, id_filename, marker_position = argv

bare = r"(^.+)(_)([a-zA-Z0-9]+)(.CEL$)"
matched = r"\3"
id_filename_dict = dd(tuple)
with open (id_filename) as id_filename_h:
    for line in id_filename_h:
        lines = line.strip().split("\t")
        #Format the filename to the same format as in the map haplotype files
        result = re.sub(bare, matched, lines[1])
        id_filename_dict[result] = lines[1] + "(" + lines[0] + ")"

marker_position_dict = dd()
with open (marker_position) as marker_position_h:
    for line in marker_position_h:
        lines = line.strip().split(",")
        marker_position_dict[lines[0]] = int(float(lines[2]))

data_haplotypes = dd(dict)
markers = list()

with open (map_file) as map_file_h:
    headers = map_file_h.readline().strip().split(",")
    for line in map_file_h:
        lines = line.strip().split(",")
        marker = lines[0]
        #All markers
        markers.append(marker)
        start = 8
        for i, (a, b) in enumerate(zip(lines[8::2], lines[9::2])):
            if a == "B":
                a = a.replace("B", "C")
            if b == "B":
                b = b.replace("B", "C")
            names = headers[start].split("-")
            true_name = names[0]
            data_haplotypes[true_name][marker] = a + " " + b
            if len(headers) >= start+2:
                start += 2
            else:
                pass
            
#Output ped file and info file.
bare = r"(^[a-zA-Z0-9]+)(_haplotypes.csv$)"
out_ped = r"\1_phased.ped"
out_info = r"\1.info"
out_ped_h = re.sub(bare, out_ped, map_file)
out_info_h = re.sub(bare, out_info, map_file)

out_ped_ho = open(out_ped_h, 'w')
out_info_ho = open(out_info_h, 'w')

#Sort marker list according to their position in the genome
selected_marker_positions = dd()
for m in markers:
    if m in marker_position_dict:
        marker_position = marker_position_dict[m]
        #print marker_position
        selected_marker_positions[marker_position] = m
#Sort the keys by number
selected_marker_positions_keys = selected_marker_positions.keys()
selected_marker_positions_keys.sort(key=int)
#Print info file
for s in selected_marker_positions_keys:
    out_info_ho.write(selected_marker_positions[s] + "\t" + str(s) + "\n")
out_info_ho.close()

#Print ped file
for sample in data_haplotypes:
    full_filename = id_filename_dict[sample]
    out_ped_ho.write("fam" + " " + full_filename + " 0 0 1 0 ")
    #print phased genotypes
    for s in selected_marker_positions_keys:
        m = selected_marker_positions[s]
        haplotype = data_haplotypes[sample][m]
        out_ped_ho.write(haplotype + " ")
    out_ped_ho.write("\n")

out_ped_ho.close()