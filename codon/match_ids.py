from sys import argv
import os, sys, re

#Takes two files and only prints out lines in the second file
#which match somewhere to the first column (ids) from the first file

script, ff, sf = argv

bare = r"(\w+)(.\w+$)"
e1 = r"\1"
e2 = e1
ff_o = re.sub(bare, e1, os.path.basename(ff))
sf_o = re.sub(bare, e2, os.path.basename(sf))
o = ff_o + "_" + sf_o
print(o)
out = open(o, 'w')

ids = {}

ff_h = open(ff)
for line in ff_h:
    fields = line.split("\t")
    gene_id = fields[0].strip()
    ids[gene_id] = 1
ff_h.close

sf_h = open(sf)
for line in sf_h:
    fields = line.split("\t")
    (match, rest) = re.split('\W+', fields[0], 1)
    if match in ids:
        out.write(line)
sf_h.close
out.close
