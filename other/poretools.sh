#!/bin/bash

#minION run stats with poretools, and conversion from FAST5 to FASTQ and FASTA.
#Input: Directory with nanopore reads in FAST5 format.

fast5=$1

poretools=/home/sobczm/bin/poretools/poretools

#Extract sequences in FASTQ format from a set of FAST5 files.
python $poretools/poretools fastq $fast5 >"${fast5}.fastq"
#Extract sequences in FASTA format from a set of FAST5 files.
python $poretools/poretools fasta $fast5 >"${fast5}.fasta"

#Get read size stats for a set of FAST5 files
python $poretools/poretools stats $fast5 >"${fast5}.stats"

#Extract the lengths and name/seq/quals from a set of FAST5 files in TAB delimited format
python $poretools/poretools tabular $fast5 >"${fast5}.tab"

#Get the qual score composition of a set of FAST5 files
python $poretools/poretools qualdist $fast5 >"${fast5}.qual"
