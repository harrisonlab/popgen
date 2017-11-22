#! /usr/bin/env python
from sys import argv
import os, sys, re
from collections import defaultdict as dd

script, annotation1, annotation2 = argv

GO_terms_all = dd(list)

with open (annotation1) as annotation1_h:
    for line in annotation1_h:
        fields = line.strip().split("\t")
        GO_terms_all[fields[0]].append(fields[3])


with open (annotation2) as annotation2_h:
    for line in annotation2_h:
        fields = line.strip().split("\t")
        GO_terms_all[fields[0]].append(fields[3])


for keys in GO_terms_all:
    all_gos = set(GO_terms_all[keys])
    all_gos = list(all_gos)
    print keys, "\t", ", ".join(all_gos)
