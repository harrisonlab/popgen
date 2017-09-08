#!/usr/bin/python

#Convert from .alt files to .csv at output by conv_cl2affy3.py script with genotypes encoded 
#as Affycodes. 

import sys
altfile = sys.argv[1]

#affy codes: 0=aa, 1=ab, 2 =bb, -1=missing
def convert_to_affy(a, b):
    if a != b:
        return '1'
    elif (a == "A" and b == "A"):
        return '0'
    elif (a == "B" and b == "B"):
        return '2'
    else:
        return '-1'

f = open(altfile)
for line in f:
    tok = line.strip().split(",")
    genotypes = list() 
    snp_id = tok[0]
    genotypes.append(snp_id)
    #Iterate over progeny genotypes:
    length_l = len(tok)
    for index in range(8,length_l,2):
        genotype = convert_to_affy(tok[index], tok[index+1])       
        genotypes.append(genotype)
    #Add the rg and ha genotypes to the end of the list.abs
    genotypes.append(convert_to_affy(tok[4], tok[5]))
    genotypes.append(convert_to_affy(tok[6], tok[7]))
    print ",".join(genotypes)