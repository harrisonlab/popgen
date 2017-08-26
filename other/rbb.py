#!/usr/bin/python
import sys
from sys import argv
#A simple script for Reciprocal Best Blast. Just what it say on the tin.abs
#Input: Species 1 vs Species 2 BLAST results, Species 2 vs Species 1 BLAST results
#Output: a list of species 1 genes in the first column with the list of species 2 genes with the matching RBB hit in the second column. 

script, s1_vs_s2, s2_vs_s1 = argv

#Blast outputs the results so that the first hit for a given query is the one with the 
#the highest score (top_candidate)

def parse_blast(my_input):
#A list to store query ids which were already seen (since we are only interested in the first oneW)
    seen = dict()
    my_output = dict()
    with open(my_input) as fh_s1:
        for line in fh_s1:
            fields = line.split()
            if fields[0] not in seen:
                #A dictionary to hold query and matching top hit.
                my_output[fields[0]] = fields[1]
                seen[fields[0]] = 1
        return my_output

s1s2 = parse_blast(s1_vs_s2)
s2s1 = parse_blast(s2_vs_s1)
#Check for best blast hits.
for k in s1s2.keys():
    hit = s1s2[k]
    if (hit in s2s1 and s2s1[hit] == k):
        print k, hit