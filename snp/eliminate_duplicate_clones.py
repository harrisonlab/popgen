#!/usr/bin/python
import sys
import re
from sys import argv
from collections import defaultdict as dd

#The script takes in a file output by the script and outputs a thinned sample id list with only one sample id per clone (cultivar) - select sample with the highest number of genotypes. 

script, input_list = argv

