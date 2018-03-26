#! /usr/bin/env python
import os, sys, re
from sys import argv
from collections import defaultdict as dd

script, one_d2, one_d = argv

keep_reads = dd(list)
with open (one_d2) as one_d2_h:
    for line in one_d2_h:
        lines = line.strip().split()
        my_id = " ".join(lines[1:])
        sequence = next(one_d2_h).strip()
        plus_sign = next(one_d2_h).strip()
        scores = next(one_d2_h).strip()
        my_true_id = line.strip()
        keep_reads[my_id] = [sequence, plus_sign, scores, my_true_id]

seen = dd()
with open (one_d) as one_d_h:
    for line in one_d_h:
        lines = line.strip().split()
        my_id = " ".join(lines[1:])
        my_true_id = line.strip()
        sequence = next(one_d_h).strip()
        plus_sign = next(one_d_h).strip()
        scores = next(one_d_h).strip()
        if my_id in keep_reads:
            print keep_reads[my_id][3]
            print keep_reads[my_id][0]
            print keep_reads[my_id][1]
            print keep_reads[my_id][2]
            seen[my_id] = 1
        else:
            print my_true_id
            print sequence
            print plus_sign
            print scores

all_1d2 = keep_reads.keys()

for my_id in all_1d2:
    if my_id in seen:
        pass
    else:
        print keep_reads[my_id][3]
        print keep_reads[my_id][0]
        print keep_reads[my_id][1]
        print keep_reads[my_id][2]