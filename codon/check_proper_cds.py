#! /usr/bin/env python
# Checks for proper start codon at the start, complete stop codon at the end
# and no in-frame stop codons
# If pass, print to the pass file, if fail, print to the fail file

from sys import argv
import os, sys, re

from Bio import SeqIO

script, fasta_file = argv

bare = r"(\w+)(.fasta)"
failed = r"\1_fail.fasta"
passed = r"\1_pass.fasta"
f = re.sub(bare, failed, fasta_file)
p = re.sub(bare, passed, fasta_file)
out_f = open(f, 'w')
out_p = open(p, 'w')

stop_codons = ["TGA", "TAG", "TAA"]

for seq_record in SeqIO.parse(fasta_file, "fasta"):
    l = len(seq_record.seq)
    #is the seq length divisable by 3 without remainder?
    r = l % 3
    #print(r)
    if r != 0:
        out_f.write(">" + seq_record.id + "\n" + str(seq_record.seq).upper() + "\n")
    else:
        start = str(seq_record.seq[0:3])
        end = str(seq_record.seq[-3:])
        t = l - 3
        is_stop_codon = lambda x : any(x[i:i+3] in stop_codons for i in range(0,t,3))
        isp = str(is_stop_codon(str(seq_record.seq)))
        is_stop_codon_r = lambda x : any(x[i:i+3] in stop_codons for i in range(0,l,3))
        ia = str(is_stop_codon_r(str(seq_record.seq)))
        #Check if sequence starts with ATG
        if start != "ATG":
            out_f.write(">" + seq_record.id + "\n" + str(seq_record.seq).upper() + "\n")
        #Check if any premature stop codon
        #Remember that elif comes straight after if block, otherwise won't compile!!
        elif isp == "True":
            out_f.write(">" + seq_record.id + "\n" + str(seq_record.seq).upper() + "\n")
        #Check if proper stop codon at the end
        elif ia == "True":
            out_p.write(">" + seq_record.id + "\n" + str(seq_record.seq).upper() + "\n")

out_f.close()
out_p.close()
