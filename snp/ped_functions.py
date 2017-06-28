from collections import defaultdict as dd
import re

def read_in_marker_definition (infile):
    count = 0
    dict_file = dd()
    with open (infile) as infile_h:
        for line in infile_h:
            lines = line.strip().split()
            dict_file[count] = lines[0]
            count += 1
    return dict_file

def read_in_marker_definition_2 (infile):
    count = 0
    dict_file = dd()
    with open (infile) as infile_h:
        for line in infile_h:
            lines = line.strip().split()
            dict_file[count] = lines[1]
            count += 1
    return dict_file

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
                    haplotype[sample_name][info_file_pos[count]] = haplo
                count += 1
    return haplotype
