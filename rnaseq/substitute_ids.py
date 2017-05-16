#!/usr/bin/python
#Read in BLAST results and substitute the hit IDs for Query IDs in a table (first column). 

from sys import argv
import os, sys, re

script, blast, table = argv

bare = r"(\w+)(.[a-zA-Z]+)"
one = r"\1_sub_ids.txt"
o = re.sub(bare, one, table)
out = open(o, 'w')

assignment = dict()

blast_h = open(blast)
for line in blast_h:
    lines = line.split()
    id = lines[1].split(".") 
    assignment[lines[0]] = id[0]
blast_h.close()

table_h = open(table)
for line in table_h:
    if lines[0] in assignment:
    	lines[0] = assignment[lines[0]]
    myString = "\t".join(lines)
    out.write(myString + "\n")
table_h.close()
out.close()