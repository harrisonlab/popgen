#!/bin/bash

#Using trimmed FASTA alignments and nucleotide substitution models identified with PartitionFinder:
#create an XML input file using BEAUTi, with StarBeast template.
#Prepare a 30 loci dataset, in addition to a 5 loci subset to compare convergence.
#Run after qlogin into a worker node (BEAST does not find BEAGLE libraries when using qsub - 
#as the BEAST package is quite fiddly, may troubleshoot it later when necessary.
#StarBeast settings:

#Substitution rate: default HKY
#Strict clock
#Species Tree Population Size: Linear with constant root
#Yule prior on species tree
#Chain length: 300 million
#Store every: 10000


path=/home/sobczm/popgen/phylogenetics/beast_runs
beast=/home/sobczm/bin/beast/BEASTv2.4.2/bin/beast
input=/<input_xml_file_here>/

cd $path
$beast -threads -1 $input
burnin=10 #percentage of states to be considered as burnin

TreeAnnotator=/home/sobczm/bin/beast/BEASTv2.4.2/bin/treeannotator

#After the run, check convergence with Tracer, summarise with TreeAnnotator and visualise with FigTree
for t in *.trees
do
output2="${t%.trees}_summary.tree"
$TreeAnnotator -heights median -burnin $burnin $t $output2
done
