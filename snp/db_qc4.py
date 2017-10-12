#!/usr/bin/python
import sys
import re
from sys import argv
from collections import defaultdict

#For cases, where istraw90 and istraw35 both contain one unique Affx, check if they match - if they don't print them out. Also, check the total number of matches.

script, alias = argv

istraw35_markers = {}
istraw90_markers = {}

with open (alias) as alias_h:
    headers = alias_h.readline().strip().split("\t")
    for line in alias_h:
        lines = line.strip().split("\t")
        affx = lines[3]
        marker = lines[4]
        if lines[1] == "istraw90":
            if marker in istraw90_markers:
                istraw90_markers[marker] = "90"
            else:
                istraw90_markers[marker] = affx
        elif lines[1] == "istraw35":
            if marker in istraw35_markers:
                istraw35_markers[marker] = "35"
            else:
                istraw35_markers[marker] = affx

counter = 0
for marker in istraw90_markers:
    if marker in istraw35_markers:
        if istraw90_markers[marker] == istraw35_markers[marker]:
            #counter += 1
            print marker, "\t", istraw35_markers[marker], "\t", istraw90_markers[marker]
        else:
            if (istraw90_markers[marker] != "90" and istraw35_markers[marker] != "35"):
                print marker, "\t", istraw35_markers[marker], "\t", istraw90_markers[marker]

#print "The total number of matching unique Affx is: ", counter