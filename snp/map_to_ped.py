#!/usr/bin/python
import sys
import re
from sys import argv
from collections import defaultdict as dd

#First argument: map based haplotypes 
#Second argument: sample_id and filename obtained from the sample table.
#Third argument: mapping of marker names to positions on the chromosome. 

script, map_file, id_filename, marker_position = argv

#Match up only the first part of the CEL file name with sample ID.
bare = r"(^.+)(_)([a-zA-Z0-9]+)(.CEL$)"
matched = r"\3"
#Dictionary to hold the sample ID.
id_filename_dict = dd(tuple)
with open (id_filename) as id_filename_h:
    for line in id_filename_h:
        lines = line.strip().split("\t")
        #Format the filename to the same format as in the map haplotype files and use as key
        result = re.sub(bare, matched, lines[1])
        #Value of the dictionary is the second field in the future ped file = CEL filename (sample_id)
        id_filename_dict[result] = lines[1] + "(" + lines[0] + ")"
#Dictionary to hold marker position on the chromosome (key: marker_id)
marker_position_dict = dd()
with open (marker_position) as marker_position_h:
    for line in marker_position_h:
        lines = line.strip().split(",")
        marker_position_dict[lines[0]] = int(float(lines[2])) #Had to do this, as due to some processing error marker positions not always integers

#dictionary of dictonaries to hold the map-based haplotypes
data_haplotypes = dd(dict)
#List of all the marker names featured in the file with map-based haplotypes.
markers = list()
#Open the custom file with map-based haplotypes
with open (map_file) as map_file_h:
    #Capture the header with sample names
    headers = map_file_h.readline().strip().split(",")
    for line in map_file_h:
        lines = line.strip().split(",")
        #Capture marker_id
        marker = lines[0]
        markers.append(marker)
        #Start with the first field featuring sample haplotypes
        start = 8
        #Iterate over maternal and paternal haplotype at the same time
        for i, (a, b) in enumerate(zip(lines[8::2], lines[9::2])):
            #Replace Bs with Cs to match ped file requirements
            if a == "B":
                a = a.replace("B", "C")
            #elif a == "A":
             #   a = a.replace("A", "C")
            if b == "B":
                b = b.replace("B", "C")
            #elif b == "A":
             #   b = b.replace("A", "C")

            #Obtain the cultivar names and sample ids so that they match up with the keys from id_filename_dict. 
            #Use as the key to data_haplotypes dictionary
            names = headers[start].split("-")
            true_name = names[0]
            #Save the haplotype for a given marker under the dictionary 
            data_haplotypes[true_name][marker] = a + " " + b
            #Move down the line onto the next sample
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

#A dictionary containing chromosomal positions of only the markers in the file with map-based haplotypes
selected_marker_positions = dd()
#Go over the markers featured in the file with map-based haplotypes
for m in markers:
    #Check that we have their position on the chromosome
    if m in marker_position_dict:
        #Retrieve the position of the marker
        marker_position = marker_position_dict[m]
        #Save the position of the marker as key and marker name as value.
        selected_marker_positions[marker_position] = m
#Save the keys to a list and sort that list by number (by value, ascending) 
selected_marker_positions_keys = selected_marker_positions.keys()
selected_marker_positions_keys.sort(key=int)
#Print info file
for s in selected_marker_positions_keys:
    out_info_ho.write(selected_marker_positions[s] + "\t" + str(s) + "\n")
out_info_ho.close()

#Print ped file
#Iterate over individuals
for sample in data_haplotypes:
    full_filename = id_filename_dict[sample]
    out_ped_ho.write("fam" + " " + full_filename + " 0 0 1 0 ")
    #print phased genotypes
    #Go over variant positions in the ascending order, as in info file
    for s in selected_marker_positions_keys:
        #Obtain the marker id to be printed
        m = selected_marker_positions[s]
        #Obtain the haplotype for a given sample and marker id
        haplotype = data_haplotypes[sample][m]
        out_ped_ho.write(haplotype + " ")
    out_ped_ho.write("\n")

out_ped_ho.close()