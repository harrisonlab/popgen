#!/usr/bin/python
import sys
import re
from sys import argv
from collections import defaultdict
#For istraw35/90 multiforms separately, need to check that no duplicates linking the same marker Affx id to the same alleles with the same flanking sequence (but MAYBE use different probes) and the same marker id - highlight them, and choose one with more genotypes available to proceed in each chip. 
script, alias, marker = argv
##Combinations which are equivalent and need to be highlighted:
#A) ref,alt,up,down - agree
#B) up, down agree; ref, alt switched
#C) ref, alt agree; up, down switched
#D) up, down switched; ref, alt switched
#E) up - reverse complement of down, down - reverse complement of up, alt - reverse complement of alt, ref - reverse complement of ref
#F) up - reverse complement of down, down - reverse complement of up, alt - reverse complement of ref, ref - reverse complement of alt

A_fh = open("A_matching.txt", 'w')
B_fh = open("B_matching.txt", 'w')
C_fh = open("C_matching.txt", 'w')
D_fh = open("D_matching.txt", 'w')
E_fh = open("E_matching.txt", 'w')
F_fh = open("F_matching.txt", 'w')

istraw35_markers = defaultdict()
istraw35_marker_ids = dict()
istraw90_markers = defaultdict()
istraw90_marker_ids = dict()

markers_dict = defaultdict()
#Read in the marker alleles, flanking data etc.
with open (marker) as marker_h:
    #headers = marker_h.readline().strip().split("\t")
    for line in marker_h:
        lines = line.strip().split("\t")
        markers_dict[lines[0]] = [lines[2:6]]

#Read in the lists of multiform markers for istraw35 and istraw90, independently.
with open (alias) as alias_h:
    #headers = alias_h.readline().strip().split("\t")
    for line in alias_h:
        lines = line.strip().split("\t")
        affx = lines[3]
        if lines[1] == "istraw90":
            istraw90_marker_ids[lines[4]] = 1
            if affx in istraw90_markers:
                istraw90_markers[affx].append(lines[4])
            else:
                istraw90_markers[affx] = [lines[4]]
        elif lines[1] == "istraw35":
            istraw35_marker_ids[lines[4]] = 1
            if affx in istraw35_markers:
                istraw35_markers[affx].append(lines[4])
            else:
                istraw35_markers[affx] = [lines[4]]

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])

#A) ref,alt,up,down - agree
def test_A(m1_ref, m2_ref, m1_alt, m2_alt, m1_up, m2_up, m1_down, m2_down, k1, k2, affx):
    if (m1_up == m2_up and m1_down == m2_down):
        if (m1_ref == m2_ref and m1_alt == m2_alt):
            A_fh.write(affx + "\t" + k1 + "\t" + k2 +  "\n")
#B) up, down agree; ref, alt switched        
def test_B(m1_ref, m2_ref, m1_alt, m2_alt, m1_up, m2_up, m1_down, m2_down, k1, k2, affx):
    if (m1_up == m2_up and m1_down == m2_down):
        if (m1_ref == m2_alt and m1_alt == m2_ref):
            B_fh.write(affx + "\t" + k1 + "\t" + k2 +  "\n")

#C) ref, alt agree; up, down switched
def test_C(m1_ref, m2_ref, m1_alt, m2_alt, m1_up, m2_up, m1_down, m2_down, k1, k2, affx):
    if (m1_ref == m2_ref and m1_alt == m2_alt):
        if (m1_up == m2_down and m1_down == m2_up):
            C_fh.write(affx + "\t" + k1 + "\t" + k2 +  "\n")


#D) up, down switched; ref, alt switched
def test_D(m1_ref, m2_ref, m1_alt, m2_alt, m1_up, m2_up, m1_down, m2_down, k1, k2, affx):
    if (m1_ref == m2_alt and m1_alt == m2_ref):
        if (m1_up == m2_down and m1_down == m2_up):
            D_fh.write(affx + "\t" + k1 + "\t" + k2 +  "\n")

#E) up - reverse complement of down, down - reverse complement of up, alt - reverse complement of alt, ref - reverse complement of ref
def test_E(m1_ref, m2_ref, m1_alt, m2_alt, m1_up, m2_up, m1_down, m2_down, k1, k2, affx):
    if (m1_ref == reverse_complement(m2_ref) and m1_alt == reverse_complement(m2_alt)):
        if (m1_up == reverse_complement(m2_down) and m1_down == reverse_complement(m2_up)):
            E_fh.write(affx + "\t" + k1 + "\t" + k2 +  "\n")

#F) up - reverse complement of down, down - reverse complement of up, alt - reverse complement of ref, ref - reverse complement of alt
def test_F(m1_ref, m2_ref, m1_alt, m2_alt, m1_up, m2_up, m1_down, m2_down, k1, k2, affx):
    if (m1_ref == reverse_complement(m2_alt) and m1_alt == reverse_complement(m2_ref)):
        if (m1_up == reverse_complement(m2_down) and m1_down == reverse_complement(m2_up)):
            F_fh.write(affx + "\t" + k1 + "\t" + k2 +  "\n")

def markers_parse(istraw90_markers, markers_dict):
    for affx in istraw90_markers:
        #Focus only on multiform markers
        if len(istraw90_markers[affx]) > 1:
            temp_dict = defaultdict()
            markers = istraw90_markers[affx]
            for m in markers:
                temp_dict[m] = markers_dict[m]
        #Compare all values against all, including reverse complement.
        #See combinations tested, above.
            for k1 in temp_dict.keys():
                for k2 in temp_dict.keys():
                        m1_ref = temp_dict[k1][0][0]
                        m2_ref = temp_dict[k2][0][0]
                        m1_alt = temp_dict[k1][0][1]
                        m2_alt = temp_dict[k2][0][1]
                        m1_up = temp_dict[k1][0][2]
                        m2_up = temp_dict[k2][0][2]
                        m1_down = temp_dict[k1][0][3]
                        m2_down = temp_dict[k2][0][3]
                        if k1 != k2:
                            test_A(m1_ref, m2_ref, m1_alt, m2_alt, m1_up, m2_up, m1_down, m2_down, k1, k2, affx)
                            test_B(m1_ref, m2_ref, m1_alt, m2_alt, m1_up, m2_up, m1_down, m2_down, k1, k2, affx)
                            test_C(m1_ref, m2_ref, m1_alt, m2_alt, m1_up, m2_up, m1_down, m2_down, k1, k2, affx)
                            test_D(m1_ref, m2_ref, m1_alt, m2_alt, m1_up, m2_up, m1_down, m2_down, k1, k2, affx)
                            test_E(m1_ref, m2_ref, m1_alt, m2_alt, m1_up, m2_up, m1_down, m2_down, k1, k2, affx)
                            test_F(m1_ref, m2_ref, m1_alt, m2_alt, m1_up, m2_up, m1_down, m2_down, k1, k2, affx)


markers_parse(istraw90_markers, markers_dict)
markers_parse(istraw35_markers, markers_dict)


A_fh.close()
B_fh.close()
C_fh.close()
D_fh.close()
E_fh.close()
F_fh.close()