#!/bin/bash

#Align minION nanopore reads using marginalign to a reference assembly, and generate
#alignment report with nanoOK.

###############Input
## First argument
#Specify the path to a directory containing the "downloaded" folder from Metrichor with the subdirectories "pass" and "fail"
## Second argument
#Specify the path to the assembly reference file, which reads are to be aligned to.

#If using barcoded data, specify the -barcoding option.

path=$1
reference=$2

#The program expects the input directory structure to reflect that output by Metrichor
#       downloaded
#            pass
#            fail
# And will output its results in the folder containing the "downloaded" subfolder

#Extract FASTA files from FAST5 input (one read -> one fasta file)
nanook extract -a -f downloaded -s $path

#Extract FASTQ files from FAST5 input (one read -> one fastq file)
nanook extract -q -f downloaded -s $path

#Align the reads with marginalign. -t argument specific number of threads to be used.
nanook align -f downloaded -s $path -r $reference -aligner marginalign -t 16
nanook analyse  -f downloaded -s $path -r $reference -aligner marginalign -t 16
