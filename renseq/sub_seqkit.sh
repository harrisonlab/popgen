#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l virtual_free=100G
#$ -l h=blacklace11.blacklace

#Sort reads in a FASTQ file by ID.
reads=$1
out=${reads%.fastq}_sort.fastq

seqkit=/home/sobczm/bin/anaconda2/bin/seqkit

$seqkit sort $reads -o $out