#!/usr/bin/python
import sys
from sys import argv
from collections import defaultdict as dd

#A script to parse the output of the rbb.py script. First argument is a plain text file with 1 column listing genes of interest, and the rest of the arguments consists of a list of 2 column tables output by rbb.py
#The script will further filter each table to output only records matching a gene in the list of genes of interest (suffix MATCH). Furthermore, for each table it will output records matching a gene in the list of genes of interest present across ALL tables and NONE of the tables.

list_genes = argv[1]
rbb_tables = argv[2:]

#read in the list of genes of interest.
goilst = dict()
with open (list_genes) as lgh:
    for line in lgh:
        goilst[line.strip()] = 1

#Iterate over the individual tables from rbb.py 
all_data = dd(lambda: dd(str))
for table in rbb_tables:
    t = open(table)
    out_m = open(table.replace(".tophits", ".match"),'w')
    out_n = open(table.replace(".tophits", ".nonmatch"),'w')
    for line in t:
        lines = line.split()
        #Check if the line matches the gene of interest - if so output
        if lines[0] in goilst:
            out_m.write(line)
            all_data[table][lines[0]] = line
    seen_ef = all_data[table].keys()
    for g in goilst:
        if g not in seen_ef:
            out_n.write(g+"\n")
    out_m.close()
    out_n.close()

#Print out lines where gene of interest present across ALL tables or none of them.
target = len(rbb_tables)
print_all = []
print_none = []
for gene in goilst:
    counter_all = 0
    counter_none = 0
    for table in all_data:
        if gene in all_data[table]:
            counter_all += 1
        else:
            counter_none += 1
    if counter_all == target:
        print_all.append(gene)
    elif counter_none == target:
        print_none.append(gene)
    else:
        pass

#print print_all

def print_lines(my_genes, my_data, my_tables, my_suffix):
    for t in my_tables:
        out = open(t.replace(".tophits", my_suffix),'w')
        for g in my_genes:
            if g in my_data[t]:
                out.write(my_data[t][g])
            else:
                out.write(g + "\n")
        out.close()

print_lines(print_all, all_data, rbb_tables, ".all")
print_lines(print_none, all_data, rbb_tables, ".none")