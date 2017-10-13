#!/usr/bin/python
import sys
import re
from sys import argv
from collections import defaultdict
#Since we have few (~13k) unique (one marker_jd -> one probe -> one Affx code) in the intersect between istraw35 and istraw90 with genotypes on both chips, would like to "rescue" some from the multiform markers. In order to do this, linking between multiform markers on istraw35 and istraw90. Matching up the same alleles, and the same flanking sequence, and keeping those markers in the analysis. Only one such pair per marker id if more available - choose the one with most genotypes on both chips. 
script, alias, marker, genotype = argv