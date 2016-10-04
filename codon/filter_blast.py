#! /usr/bin/env python
from __future__ import division
import os, sys, re, argparse


ap = argparse.ArgumentParser()
ap.add_argument('--i',required=True,type=str,help='Input Blast results')
#Filter by Percentage Identity (3rd column)
#default 90%
ap.add_argument('--id',required=False,type=float,help='Minimum percentage identity',default=90)
ap.add_argument('--cov',required=False,type=float,help='Minimum percentage subject coverage', default=80)
#Filter by Subject Coverage (length of the alignment -> 4th column / subject length -> 15th column )
#default 80%
args = ap.parse_args()

blast=args.i
id_th=args.id
cov_th=args.cov

o = blast + "_filtered"

blast_out = open(o, 'w')

blast_h = open(blast)
for line in blast_h:
    fields = line.split()
    print(fields[14])
    coverage_per = (int(fields[3]) / int(fields[14])) * 100
    #Percentage Identity
    if float(fields[2]) < id_th:
        continue
    elif coverage_per < cov_th:
        continue
    else:
        blast_out.write(line)

blast_h.close
blast_out.close
