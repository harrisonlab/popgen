#!/bin/bash
input=/home/sobczm/popgen/renseq/input/transcriptomes/really_really_final_baits
scripts=/home/sobczm/bin/popgen/renseq
#Final strategy for bait design involving RLP, RLK and SRLK genes.
#Cluster the nucleotide (protein?) sequences at different levels of sequence identity and take
#represenative contig. Design the bait to the contigs at 5x density and count up the total number
#of baits. Stop when reaching >100,000 baits in total.

mkdir $input
#Collect all the results required for final bait design here following rgaugury.sh script execution.
