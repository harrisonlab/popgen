#!/bin/bash

#Using trimmed FASTA alignments and nucleotide substitution models identified with PartitionFinder:
#create an XML input file using BEAUTi, with StarBeast template.
#Prepare a 20 loci dataset, in addition to a 5 loci subset to compare convergence.
#Run after qlogin into a worker node (BEAST does not find BEAGLE libraries when using qsub -
#as the BEAST package is quite fiddly, may troubleshoot it later when necessary.

#StarBeast settings:

#Substitution rate: default GTR
#Strict clock
#Species Tree Population Size: Constant
#Yule prior on species tree

path=/home/sobczm/popgen/phylogenetics/clock/CDS_genomes/beast_runs/candidates/trimmed
beast=/home/sobczm/bin/beast/BEASTv2.4.2/bin/beast
input=/<input_xml_file_here>/

cd $path
$beast -threads -1 $input

TreeAnnotator=/home/sobczm/bin/beast/BEASTv2.4.2/bin/treeannotator

#After the run, check convergence with Tracer, summarise with TreeAnnotator and visualise with FigTree
burnin=10 #percentage of states to be considered as burnin
for t in *.trees
do
output2="${t%.trees}_summary.tree"
$TreeAnnotator -heights median -burnin $burnin $t $output2
done
