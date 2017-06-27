#!/usr/bin/python
import sys
import re
from sys import argv
from collections import defaultdict as dd

#First argument: map based haplotypes 
#Second argument: sample_id and filename obtained from the sample table.abs
#Third argument: mapping of marker names to positions on the chromosome. 

script, map_file, id_filename, marker_position = argv

id_filename_dict = dd()
with open (id_filename) as id_filename_h:
    for line in id_filename_h:
        lines = line.strip().split("\t")
        id_filename_dict[lines[0]] = lines[1]

marker_position_dict = dd()
with open (marker_position) as marker_position_h:
    for line in marker_position_h:
        lines = line.strip().split(",")
        marker_position_dict[lines[0]] = lines[2]

