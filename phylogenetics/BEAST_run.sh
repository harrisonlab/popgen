#!/bin/bash
#Using trimmed FASTA alignments and nucleotide substitution models identified with PartitionFinder:
#create an XML input file using BEAUTi, with StarBeast template.
# Beauti: /home/sobczm/bin/beast/BEASTv2.4.2/bin/beauti
#Prepare a 30 loci dataset, in addition to a 5 loci subset to compare convergence.
#Run after qlogin into a worker node (BEAST does not find BEAGLE libraries when using qsub - 
#as the BEAST package is quite fiddly, may troubleshoot it later when necessary.

#StarBeast settings used here:

#Substitution rate: default HKY
#Strict clock
#Species Tree Population Size: Linear with constant root
#Yule prior on species tree
#Chain length: 300 million (this may vary, change run convergence with Tracer during the run to establish the number of iterations required
#Tracer: /home/sobczm/bin/beast/Tracer_v1.6/bin/tracer
#some runs may never converge)
#Store every: 10000

#First login into the selected node. Start running the script below to reserve the resources on selected node (need to change the node in the file header)
#https://github.com/harrisonlab/popgen/blob/master/other/reserve_cluster.sh

path=/home/sobczm/popgen/phylogenetics/beast_runs
beast=/home/sobczm/bin/beast/BEASTv2.4.2/bin/beast
input=/<input_xml_file_here_created_by_beauti>/

cd $path
$beast -threads -1 $input

#After the run, check convergence with Tracer, summarise the final tree with TreeAnnotator 
burnin=10 #percentage of states to be considered as burnin

TreeAnnotator=/home/sobczm/bin/beast/BEASTv2.4.2/bin/treeannotator
for t in *.trees
do
output2="${t%.trees}_summary.tree"
$TreeAnnotator -heights median -burnin $burnin $t $output2
done

#Visualise and beautify the final tree (suffix "summary") with FigTree
/home/sobczm/bin/FigTree_v1.4.2/bin/figtree
