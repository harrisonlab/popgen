#!/usr/bin/python
import sys
import re
from sys import argv
from collections import defaultdict as dd
#See what crown rot samples are found somewhere in db: tables clone and sample, and which are not. 
script, crown_rot_tb, clones_tb, samples_tb = argv
#If the cultivar name is just a string with no numbers, lowercase it and compare to everything. 
#If the cultivar name contains numbers just grep for numbers.

#Ignore the id column in the tables.

crown_rot = dd()
#Read in the list of the cultivars in the crown rot GWAS table.
with open (crown_rot_tb) as crown_rot_tb_h:
    for line in crown_rot_tb_h:
        fields = line.split("\t")
        #Convert to lower-case:
        search_token = fields[0].lower()
       # if any(char.isdigit() for char in search_token):
            #strip all the letters:
            #search_token = re.sub("\D", "", search_token)
            #crown_rot[search_token] = fields[0]
        #else:
        crown_rot[search_token] = fields[0]

#read in the clones table
clones = dd()
with open (clones_tb) as clones_tb_h:
    for line in clones_tb_h:
       # lines = line.split("\t")
        #line = "\t".join(lines[1:])
        clones[line.lower()] = 1

#read in the samples table
samples = dd()
with open (samples_tb) as samples_tb_h:
    for line in samples_tb_h:
        #lines = line.split("\t")
        #line = "\t".join(lines[1:])
        samples[line.lower()] = 1

#print the results
tally = 0
nonmatched = list()
for my_id in crown_rot.keys():
    flag = 0
    for my_clone in clones.keys():
        if re.search(my_id, my_clone):
            flag = 1
            print my_id,
            print my_clone,
    for my_sample in samples.keys():
        if re.search(my_id, my_sample):
            flag = 1
            print my_id,
            print my_sample,
    if flag == 0:
        nonmatched.append(my_id)
    else:
        tally += 1

#Print non-matched samples at the end of the file:
for n in nonmatched:
    print crown_rot[n]
print tally
