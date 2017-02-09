#!/bin/bash
input=/home/sobczm/popgen/phylogenetics/beast_runs/loci30
phyutility=/home/sobczm/bin/phyutility

#First prune the consensus tree
java -jar $phyutility/phyutility.jar -pr -in $input/30loci_1471365593665_summary.tree -out $input/30loci_1471365593665_summary_pruned.tree -names 55 A1-2 D2 HB6 proliferatum

#Secondly prune all the replicate runs trees
java -jar $phyutility/phyutility.jar -pr -in $input/30loci_1471365593665.trees -out $input/30loci_1471365593665_pruned.trees -names 55 A1-2 D2 HB6 proliferatum

#Thirdly, calculate the new bootstrap supports
java -jar $phyutility/phyutility.jar -ts -in $input/30loci_1471365593665_pruned.trees -tree $input/30loci_1471365593665_summary_pruned.tree -out $input/30loci_1471365593665_summary_pruned_bs.tree 