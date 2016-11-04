#! /usr/bin/env python
from sys import argv
import os, sys, re

#Well, it turns out that despite filtering out unmapped reads and their mates
#the number of reads in forward and reverse files is still not equal.
#To combat that, run this script twice: forward versus reverse. and reverse vs forward
#reads to keep only reads present in both files at the same time!

script, one_rs, second_rs, output_reads = argv

#First two arguments: input
#Last two arguments: output

header_regex = r"(.+)(/\d)$"
set_reads = dict()

one_rs_h = open(one_rs)
for line in one_rs_h:
    if line.startswith("@"):
        result = re.search(header_regex, line)
        set_reads[result.group(1)] = 1

one_rs_h.close

second_rs_h = open(second_rs)
output_reads_h = open(output_reads, 'w')

for line in second_rs_h:
    if line.startswith("@"):
        result = re.search(header_regex, line)
        if result.group(1) in set_reads:
            output_reads_h.write(line)
            for i in range(3):
                next_l = second_rs_h.next()
                output_reads_h.write(next_l)
        else:
            pass

output_reads_h.close
second_rs_h.close
