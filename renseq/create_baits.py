#!/usr/bin/python

'''
create a set of bait sequence from the fasta file
output as new fasta file
'''

from Bio import SeqIO
from merge_overlaps_funcs import *
import argparse,sys

ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inp',required=True,nargs='+',type=str,help='input sequences FASTA file(s)')
ap.add_argument('--size',default=120,type=int,help='bait size (bps)')
ap.add_argument('--coverage',default=2,type=int,help='bait coverage')
ap.add_argument('--out',default='STDOUT',type=str,help='output sequences FASTA')
conf = ap.parse_args()

#open output file
if conf.out == 'STDOUT':
    fout = sys.stdout
else:
    fout = open(conf.out,'wb')

ct = 0

#load each fasta, output baits
for fname in conf.inp:
    for rec in SeqIO.parse(fname, "fasta"):
        uid = rec.id

        #convert to upper case
        #break at Ns
        seq_list = cleanup_seq(str(rec.seq))

        for seq in seq_list:
            #calc positional shift between each 1x coverage tiling set
            shift = int(round(conf.size/conf.coverage))

            for i in xrange(conf.coverage):
                for bait in generate_baits(seq[i*shift:],conf.size):
                    fout.write('>'+uid+'::'+str(ct) + '\n')
                    fout.write(bait + '\n')
                    ct += 1

if conf.out != 'STDOUT':
    fout.close()
