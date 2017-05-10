#! /usr/bin/env python
import os
import re
import sys
from collections import defaultdict
from sys import argv

from Bio import SeqIO

#Scan one strand for telomere repeats
###Arguments:
#First: Input file with sequence(s) to be scanned:
#Second: Expected repeat motif sequence (5'-> 3') in upper case
#Third: Which strand to scan?
#Output: .bed file listing matches: contig/start/end/strand

script, input_file, motif, strand = argv

bare = r"(\w+)(.fa$|.fasta$)"
out_sub = r"\1"
o = re.sub(bare, out_sub, input_file)
out_name = o + "_" + motif + "_" + strand + "strand" + ".bed"
out_name2 = "contig_lengths.txt"

out = open(out_name, 'w')
out2 = open(out_name2, 'w')

def scan_sequence(seq_name, sequence):
    a = range(len(sequence))
        #print(a)
    for n in a:
        end_l = n + len(motif)
        #Check if not going beyond the length of sequnece
        if end_l < len(sequence):
            #Check if motif match
            if sequence[n:end_l] == motif:
                out.write(seq_name + "\t" + str(n) + "\t" + str(end_l) + "\t" + strand + "\n")

print ("Searching for motif: " + motif + " on " + strand + " strand")
for seq_record in SeqIO.parse(input_file, "fasta"):
    if strand == "-":
        seq_revcom = seq_record.seq.reverse_complement()
        seq_upper = str(seq_revcom.upper())
    else: 
        seq_upper = str(seq_record.seq).upper()
    seq_name = str(seq_record.id) 
    con_len = len(seq_upper)
    #Print the length of the sequence to screen
    print ("Scanning " + seq_name + " whose total length is " + str(con_len))
    out2.write(seq_name + "\t" + str(con_len) + "\n")
    scan_sequence(seq_name, seq_upper)

out2.close()
out.close()