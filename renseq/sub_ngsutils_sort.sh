#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l virtual_free=100G
#$ -l h=blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace11.blacklace

reads=$1
out=${reads%.fastq}_sort.fastq

ngsutils=/home/sobczm/bin/ngsutils/bin

$ngsutils/fastqutils sort $reads >$out