#! /usr/bin/env python
from sys import argv
import os, sys, re

# First read in mapped reads from TopHat analysis and save read names
# (minus the read orientation) into memory.
# Then start reading in unmapped reads from TopHat analysis. If the read name
# (minus the read orientation) has not been previously seen, print it to file
# with forward or reverse reads. Afterwards sort the reads in bash.

script, mapped, unmapped, forward_reads, reverse_reads = argv

#First two arguments: input
#Last two arguments: output

header_regex = r"(.+)(/\d)$"
unmapped_reads = dict()

mapped_h = open(mapped)
for line in mapped_h:
    if line.startswith("@"):
        result = re.search(header_regex, line)
        unmapped_reads[result.group(1)] = 1

mapped_h.close

unmapped_h = open(unmapped)
forward_h = open(forward_reads, 'w')
reverse_h = open(reverse_reads, 'w')

for line in unmapped_h:
    if line.startswith("@"):
        result = re.search(header_regex, line)
        if result.group(1) in unmapped_reads:
            pass
        elif result.group(2) == '/1':
            forward_h.write(line)
            for i in range(3):
                next_l = unmapped_h.next()
                forward_h.write(next_l)
        elif result.group(2) == '/2':
            reverse_h.write(line)
            for i in range(3):
                next_l = unmapped_h.next()
                reverse_h.write(next_l)
        else:
            pass

unmapped_h.close()        
reverse_h.close()
