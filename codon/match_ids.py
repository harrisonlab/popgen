from sys import argv
import os, sys, re

#Takes two files and only prints out lines in the second file
#which match somewhere to the first column (ids) from the first file
#The second file contains the lines from the first file not having an id match.

script, ff, sf = argv

bare = r"(\w+)(.\w+$)"
e1 = r"\1"
e2 = e1
ff_o = re.sub(bare, e1, os.path.basename(ff))
sf_o = re.sub(bare, e2, os.path.basename(sf))
o = ff_o + "_" + sf_o
o2 = ff_o + "_" + sf_o + "_remainder"
out = open(o, 'w')
out2 = open(o2, 'w')

ids = {}

ff_h = open(ff)
for line in ff_h:
    fields = line.split()
    gene_id = fields[0].strip()
    match = re.split('\W+', gene_id, 1)
    #temp_id = "Fus2_" + match[0]
    ids[match[0]] = 1
ff_h.close

sf_h = open(sf)
for line in sf_h:
    fields = line.split()
    match = re.split('\W+', fields[0], 1)
    if match[0] in ids:
        out.write(line)
    else:
        out2.write(line)
sf_h.close
out.close
