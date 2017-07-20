#!/usr/bin/python

#
# convert raw affy calls, one marker per line, into PED format
#

import sys
import numpy as np

inpfile = sys.argv[1]

#input is one row per marker as: markerid segtype phase genotype...
data = np.transpose(np.genfromtxt(inpfile,dtype=object,delimiter='\t'))
#after transpose one row per individual
#first row is markerid... etc

#convert segtype line into two parental genotype call lines
#for maternal and paternal parents
conv =\
{
    '-1':'0 0',
    '0':'A A',
    '1':'A B',
    '2':'B B',
}

ct = 1

for row in data[1:]:
    print 'fam %s 0 0 1 0 '%row[0] + ' '.join([conv[x.strip()] for x in row[1:]])
    ct += 1
