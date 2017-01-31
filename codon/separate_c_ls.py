#! /usr/bin/env python
from sys import argv
import os, sys, re

#Split the input annotation table into genes on the core and lineage-specific cds_to_chromosome_coords

script, input_table = argv

bare = r"(\w+)(.txt$)"
e1 = r"\1"
ff_o = re.sub(bare, e1, input_table)
o_ls = ff_o + "_ls.txt"
o_c = ff_o + "_c.txt"
o_ps = ff_o + "_ps.txt"
o_ss = ff_o + "_ss.txt"
out = open(o_ls, 'w')
out2 = open(o_c, 'w')
out3 = open(o_ps, 'w')
out4 = open(o_ss, 'w')

lineage_specific = ["contig_10_pilon", "contig_14_pilon", "contig_16_pilon", "contig_18_pilon",
"contig_19_pilon", "contig_20_pilon", "contig_21_pilon", "contig_22_pilon"]
species_specific = ["contig_14_pilon", "contig_18_pilon", "contig_20_pilon", "contig_22_pilon"]
pathogen_specific = ["contig_10_pilon", "contig_16_pilon",  "contig_19_pilon", "contig_21_pilon"]
core = ["contig_1_pilon", "contig_2_pilon", "contig_3_pilon", "contig_4_pilon",
"contig_5_pilon", "contig_6_pilon", "contig_7_pilon", "contig_8_pilon", "contig_9_pilon",
"contig_11_pilon", "contig_12_pilon", "contig_13_pilon",
"contig_15_pilon", "contig_17_pilon"]

ff_h = open(input_table)
for line in ff_h:
    fields = line.split()
    contig_id = fields[1].strip()
    if contig_id in lineage_specific:
        out.write (line)
        if contig_id in species_specific:
            out4.write (line)
        elif contig_id in pathogen_specific:
            out3.write (line)
    elif contig_id in core:
        out2.write (line)

out.close()
out2.close()
ff_h.close()
