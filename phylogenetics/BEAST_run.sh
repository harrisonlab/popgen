#!/bin/bash

#Using trimmed FASTA alignments and nucleotide substitution models identified with PartitionFinder:
#create an XML input file using BEAUTi, with StarBeast template.
#Prepare a 30 loci dataset, in addition to a 5 loci subset to compare convergence.
#Run after qlogin into a worker node (BEAST does not find BEAGLE libraries when using qsub - 
#as the BEAST package is quite fiddly, may troubleshoot it later when necessary.


path=/home/sobczm/popgen/phylogenetics/beast_runs
beast=/home/sobczm/bin/beast/BEASTv2.4.2/bin/beast
input=/<input_xml_here>/

cd $path
$beast -threads -1 $input

#After the run, check convergence with Tracer, visualise the trees with DensiTree and TreeAnnotator-Figtree.
