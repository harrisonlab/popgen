#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 10
#$ -l virtual_free=10G
#$ -l h=blacklace11.blacklace

#Sort reads in a FASTQ file by ID.
reads=$1
out=${reads%.fastq}_sort.fastq

fastqt=/home/sobczm/bin/fastq-tools-0.8/bin

$fastqt/fastq-sort --id $reads > $out