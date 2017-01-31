from sys import argv
import os, sys, re
from collections import defaultdict

#Parse the summary table output by detect_duplications.py to find
#out which groups of contigs the queries and the subjects belong to:
#core, core* or lineage-specific (for contigs 1-22) in the Fus2 genome.

script, gene_table, summary = argv
contigs = defaultdict(str)

lineage_specific = ["contig_10_pilon", "contig_14_pilon", "contig_16_pilon", "contig_18_pilon",
"contig_19_pilon", "contig_20_pilon", "contig_21_pilon", "contig_22_pilon"]
core = ["contig_1_pilon", "contig_2_pilon", "contig_3_pilon", "contig_4_pilon",
"contig_5_pilon", "contig_6_pilon", "contig_7_pilon", "contig_8_pilon"]
core_asterisk = ["contig_9_pilon", "contig_11_pilon", "contig_12_pilon", "contig_13_pilon",
 "contig_15_pilon", "contig_17_pilon"]

c_ca = 0
c_c = 0
c_ls = 0

ca_c = 0
ca_ca = 0
ca_ls = 0

ls_c = 0
ls_ca = 0
ls_ls = 0


ff_h = open(gene_table)
for line in ff_h:
    fields = line.split()
    contigs[fields[0]] = fields[1]
ff_h.close

sf_h = open(summary)
for line in sf_h:
    fields = line.split("\t")
    match = re.split('\W+', fields[0], 1)
    #Get all hits
    tandem = fields[2].split(",")
    segmental = fields[4].split(",")
    print(segmental[0])
    #print(segmental)
    all_hits = tandem + segmental
    print (all_hits)
    if match[0] in contigs:
        for a in all_hits:
            match_a = re.split('\W+', a, 1)
            if ((contigs[match[0]] in core) and (contigs[match_a[0]] in core_asterisk)):
                c_ca += 1
            elif ((contigs[match[0]] in core) and (contigs[match_a[0]] in core)):
                c_c += 1
            elif ((contigs[match[0]] in core) and (contigs[match_a[0]] in lineage_specific)):
                c_ls += 1
            elif ((contigs[match[0]] in core_asterisk) and (contigs[match_a[0]] in core)):
                ca_c += 1
            elif ((contigs[match[0]] in core_asterisk) and (contigs[match_a[0]] in core_asterisk)):
                ca_ca += 1
            elif ((contigs[match[0]] in core_asterisk) and (contigs[match_a[0]] in lineage_specific)):
                ca_ls += 1
            elif ((contigs[match[0]] in lineage_specific) and (contigs[match_a[0]] in core)):
                ls_c += 1
            elif ((contigs[match[0]] in lineage_specific) and (contigs[match_a[0]] in core_asterisk)):
                ls_ca += 1
            elif ((contigs[match[0]] in lineage_specific) and (contigs[match_a[0]] in lineage_specific)):
                ls_ls += 1

sf_h.close

print ("Query core & Subject core " + str(c_c) + "\n")
print ("Query core & Subject core asterisk " + str(c_ca) + "\n")
print ("Query core & Subject lineage-specific " + str(c_ls) + "\n")

print ("Query core asterisk & Subject core " + str(ca_c) + "\n")
print ("Query core asterisk & Subject core asterisk " + str(ca_ca) + "\n")
print ("Query core asterisk & Subject lineage-specific " + str(ca_ls) + "\n")

print ("Query lineage_specific & Subject core " + str(ls_c) + "\n")
print ("Query lineage_specific & Subject core asterisk " + str(ls_ca) + "\n")
print ("Query lineage_specific & Subject lineage-specific " + str(ls_ls) + "\n")
