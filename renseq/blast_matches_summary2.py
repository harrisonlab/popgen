#!/usr/bin/python
import sys
from sys import argv
from collections import defaultdict as dd

length_genes = argv[1]
length_rgenes = argv[2]
rbb_tables = argv[3]

#read in the length of subject genes of interest.
goilst = {}
with open (length_genes) as lgh:
    for line in lgh:
        lines = line.split(",")
        goilst[lines[0]] = lines[1]

#read in the length of our assembled query genes of interest.
goilstr = {}
with open (length_rgenes) as lghr:
    for line in lghr:
        lines = line.split(",")
        goilstr[lines[0]] = lines[1]

names=["qseqid", "sseqid", "pident", "aln_length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "sstrand"]

def generate_range(start, end):
    if int(end) > int(start):
        return set(range(int(start), int(end)+1))
    else:
        return set(range(int(end), int(start)+1))

three_utrs = dict()
five_utrs = dict()

def define_utr(fields):
    if (int(fields[8]) == 1 or int(fields[9]) == 1):
        futrl = int(goilstr[fields[0]]) - int(fields[6])
        five_utrs[fields[0]] = futrl
    elif (int(fields[8]) >= int(goilst[fields[1]]) - 3 or int(fields[9]) >= int(goilst[fields[1]]) - 3):
        rutrl = futrl = int(goilstr[fields[0]]) - int(fields[7])
        three_utrs[fields[0]] = rutrl


seen_query = dict()
query_range = dd(set)
hit_range = dd(set)
with open(rbb_tables) as fh_s1:
        for line in fh_s1:
            fields = line.split()
            if fields[0] in seen_query:
                if fields[1] == seen_query[fields[0]]:
                    queryr = generate_range(fields[6], fields[7])
                    subr = generate_range(fields[8], fields[9])
                    diff_q = queryr.difference(query_range[fields[0]])
                    diff_s = subr.difference(hit_range[fields[0]])
                    query_range[fields[0]].update(diff_q)
                    hit_range[fields[0]].update(diff_s)
                    #Define UTRs
                    define_utr(fields)
                else:
                    pass
            else:
                seen_query[fields[0]] = fields[1]
                queryr = generate_range(fields[6], fields[7])
                query_range[fields[0]].update(queryr)
                subr = generate_range(fields[8], fields[9])
                hit_range[fields[0]].update(subr)
                #Define UTRs
                define_utr(fields)
                
#Sum the best HSPs
#Take whatever is smallest - query or hit coverage. 
for query in seen_query:
    print query, seen_query[query],
    if hit_range[query] > query_range[query]:
        length = len(query_range[query])
    else:
        length = len(hit_range[query])
    #Print out the query name, the hit name, alignment length, query gene length, subject gene length and % coverage
    print length, goilstr[query].strip(), goilst[seen_query[query]].strip(), float(length/float(goilst[seen_query[query]]))

#Print the UTR results: query name, hit name, UTR length, query gene length, subject gene length
def print_utr(my_results, output_file):
    out = open(output_file, 'w')
    for query in my_results:
        out.write(query + "\t" + seen_query[query] + "\t" + str(my_results[query]) + "\t" + str(goilstr[query].strip()) + "\t" + str(goilst[seen_query[query]]))
    out.close()

out5 = rbb_tables + "_5utr"
out3 = rbb_tables + "_3utr"
print_utr(five_utrs, out5)
print_utr(three_utrs, out3)