import sys
import re
from sys import argv
from collections import OrderedDict as odd
#Add sample id to the master strawberry table. Only perfect case-blind matches accepted. 
script, old_table, new_table = argv

crown_rot = dict()
#Read in the list of the cultivars in the crown rot GWAS table.
with open (old_table) as old_table_h:
    for line in old_table_h:
        fields = line.strip().split("\t")
        #Convert to lower-case:
        if len(fields) > 1:
            crown_rot[fields[0]] = fields[1]
        else:
            crown_rot[fields[0]] = ""

error = open("error.log", 'w')
with open (new_table) as new_table_h:
    for line in new_table_h:
        fields = line.strip().split("\t")
        if fields[0] in crown_rot:
            print fields[0], "\t", crown_rot[fields[0]]
        else:
            print fields[0]
            error.write(fields[0] + "\n")
error.close()