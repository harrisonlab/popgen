#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l virtual_free=15G
#$ -l h=blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

#Usage: qsub sub_pairfq.sh forward_reads reverse_reads

#The script accepts FASTQ files: forward and reverse reads.
#It re-pairs paired-end reads (output suffix: _f_p.fastq and _r_p.fastq) which are out of sync, and also outputs two files with singletons (output suffix: _f_s.fastq and _r_s.fastq), which are subsequently gunzipped (compressed).

forward=$1
reverse=$2

pairfq=/home/sobczm/bin/Pairfq/bin/pairfq

perl $pairfq makepairs -c gzip -f $forward -r $reverse -fp ${forward%.fastq}_p.fastq -rp ${reverse%.fastq}_p.fastq -fs ${forward%.fastq}_s.fastq -rs ${reverse%.fastq}_s.fastq
