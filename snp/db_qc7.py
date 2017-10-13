#!/usr/bin/python
import sys
import re
from sys import argv
from collections import defaultdict
#In addition, for istraw90 only dataset, look at multiform markers ids present only on that chip, and pick one with the most genotypes.
script, alias, marker, sample = argv