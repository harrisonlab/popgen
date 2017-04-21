#! /usr/bin/env python
from sys import argv
import os, sys, re
from collections import defaultdict
from Bio import SeqIO

###Arguments:
#First: Input file with sequence(s) to be scanned:
#Second: 5' motif flank:
#Third: 3' motif flank:
#Fourth: What is the gap length range (e.g. 0,20):
#Fifth: motif name

script, input_file, flank_left, flank_right, gap_range, motif_name = argv

gaps = gap_range.split(",")
gaps_range = range(int(gaps[0]), int(gaps[1])+1)
l_flan_len = len(flank_left)
r_flan_len = len(flank_right)
motifs = defaultdict(lambda: defaultdict(dict))


def scan_sequence(seq_name, sequence, strand):
    a = range(len(sequence))
        #print(a)
    for n in a:
        end_l = n + l_flan_len
        #Check if not going beyond the length of sequnece
        if end_l < len(sequence):
            #Check if motif match
            if sequence[n:end_l] == flank_left:
              #Look for the right part of the motif
               for b in gaps_range:
                   end_r = end_l + b + r_flan_len
                   #Check if not moving beyond the end of sequence
                   if end_r < len(sequence):
                       #Check if motif_match
                       start_l = end_l + b
                       if sequence[start_l:end_r] == flank_right:
                           #Save the motif
                           ids = str(n) + "_" + str(end_r)
                           if strand == "-":
                               end_r = len(sequence) - end_r
                               n = len(sequence) - n
                               motifs[seq_name][ids] = [end_r, n, b, strand]
                           if strand == "+":
                               motifs[seq_name][ids] = [n, end_r, b, strand]

def parse_input(input_file):
    for seq_record in SeqIO.parse(input_file, "fasta"):
        seq_upper = str(seq_record.seq).upper()
        seq_name = str(seq_record.id)
        scan_sequence(seq_name, seq_upper, "+")
        seq_revcom = seq_record.seq.reverse_complement()
        scan_sequence(seq_name, str(seq_revcom.upper()), "-")

def print_output(motif_name):
    out_name = motif_name + "_motif_list.txt"
    out_name2 =  motif_name + "_motif_summary.txt"
    out = open(out_name, 'w')
    out2 = open(out_name2, 'w')
    out.write("5' motif flank: " + flank_left + "\n")
    out.write("3' motif flank: " + flank_right + "\n")
    out.write("Gap length range: " + gap_range + "\n")
    out.write("sequence name" + "\t" +  "match start" + "\t" + "match end" + "\t" + "gap length" + "\t" + "strand" + "\n")
    for j in motifs:
        for k in motifs[j]:
             out.write(j + "\t" + str(motifs[j][k][0]+1) + "\t" + str(motifs[j][k][1]) + "\t" + str(motifs[j][k][2]) + "\t" + str(motifs[j][k][3]) + "\n")
    for j in motifs:
        motif_no = len(motifs[j])
        out2.write(j + "\t" + str(motif_no) + "\n")
    out.close()
    out2.close()


parse_input(input_file)
print_output(motif_name)