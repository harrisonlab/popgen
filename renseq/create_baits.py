#!/usr/bin/python

'''
create a set of bait sequences from the fasta file(s)
output as a new fasta file
'''

import argparse,sys

ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inp',required=True,nargs='+',type=str,help='input sequences FASTA file(s)')
ap.add_argument('--size',default=120,type=int,help='bait size (bps)')
ap.add_argument('--coverage',default=2,type=int,help='bait coverage')
ap.add_argument('--out',default='STDOUT',type=str,help='output sequences FASTA')
conf = ap.parse_args()

from Bio import SeqIO
#from merge_overlaps_funcs import *

def generate_baits(seq,size):
    '''
    chop seq up into tiled baits of the required size
    exclude the final bait if it's not full length
    '''

    nseq = len(seq)
    for i in xrange(0,nseq,size):
        bait = seq[i:i+size]
        if len(bait) != size: continue #exclude partial baits
        yield bait

def cleanup_seq(seq):
    '''
    convert to upper case
    check for unexpected symbols
    ##split at Ns
    '''

    #uppercase
    seq = seq.upper()

    #check for unexpected characters
    allow = 'ATCGN' #allowed symbols
    chk = [x for x in seq if x in allow]
    assert len(chk) == len(seq)

    #split at Ns
    #seq_list = [x for x in seq.split('N') if x != '']

    return seq_list

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
