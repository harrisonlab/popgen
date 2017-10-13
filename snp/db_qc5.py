#!/usr/bin/python
import sys
import re
from sys import argv
from collections import defaultdict
#For istraw35/90 multiforms separately, need to check that no duplicates linking the same marker id to the same alleles with the same flanking sequence but MAYBE use different probes - highlight them, and choose one with more genotypes available to proceed in each chip. 

script, alias, marker, sample = argv