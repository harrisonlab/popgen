#!/usr/bin/python

#
# convert map positions to SHAPEIT's "genetic map" file
# which specifies the recombination rate as a function of physical and genetic map posn
# we used vesca physical positions, assuming a uniform crossover rate (unrealistic)
# then scaled to make total map distance comparable to rgxha map
# => 3.125 cM per megabase

import sys

inpfile = sys.argv[1]

#convert to integer

f = open(inpfile)
#f.readline() #skip header
print 'pposition rrate gposition'
for i,line in enumerate(f):
    tok = line.strip().split()
    cm = float(tok[2])
    bp = int(tok[3])
    print bp,3.333,cm
f.close()