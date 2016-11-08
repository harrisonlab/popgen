#!/bin/bash
input=/home/sobczm/popgen/clock/DNA_genomes/promoters
scripts=/home/sobczm/bin/popgen/clock/motif_discovery
meme=/home/sobczm/bin/meme_4.11.2/bin
weeder=/home/sobczm/bin/weeder/weeder2

cd $input
#Repeat-masking of repetitive sequences already done in the downloaded genomes.

#Find known motiffs in N. crassa
#Cbox in frequency promoter
#CGAT(N)CCGCT

#N: 0-3 bp or more
# Use GLAM2Scan


#ACE element in ccg2
#The core ACE sequence binding site (Bell-Pedersen 2001), from -107,
#similar a bit to a sequence in other ccgs
#AACTTGGCCAAGTT
# Use FIMO
$meme/iupac2meme AACTTGGCCAAGTT >core_ace.txt
qsub $scripts/sub_fimo.sh core_ace.txt Neurospora_crassa.NC12.dna_rm.toplevel_promoters_1000.fasta core_ace


#The full ACE containing element (68 bp, between -118 and -50 of the cc2 TSS)
#(Bell-Pedersen 2001)
#GAATACCGGAGAACTTGGCCAAGTTTGATGGACGAAGTCTTCAAACACAGCGTTGGATTGAGGTCCAA
# Use FIMO
$meme/iupac2meme GAATACCGGAGAACTTGGCCAAGTTTGATGGACGAAGTCTTCAAACACAGCGTTGGATTGAGGTCCAA >full_ace.txt
qsub $scripts/sub_fimo.sh full_ace.txt Neurospora_crassa.NC12.dna_rm.toplevel_promoters_1000.fasta full_ace
#The top match is located in the protein EAA34064, which is gene NCU08457, which is
#ccg2, which agrees with expectations. It is on plus strand, but situated
#between 767 and 834 bp from the beginning of a 1000 bp promoter, ie.
#between -233 and -166 bp from TSS. This agrees with GFF annotation on Ensembl
#but not with the experimental data (see above). Also, the hit is not perfect but the best by an order of magnitude
