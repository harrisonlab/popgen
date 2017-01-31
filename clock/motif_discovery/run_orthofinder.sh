#!/bin/bash
#Input: directory containing all the FASTA files to be used in the analysis,
# with each containing protein sequences for one species.

#Output: Tables showing summary stats (Statistics_Overall.csv, Statistics_PerSpecies.csv, Orthogroups_SpeciesOverlaps.csv)
#and actual orthogroup assignment (Orthogroups.csv, Orthogroups.txt)
#ALSO: gene trees for each orthogroup in the Newick format (folder Trees)

dir=$1

orthofinder=/home/sobczm/bin/OrthoFinder-1.0.7/orthofinder
anaconda=/home/sobczm/bin/anaconda2/bin/python

#16 threads used
$anaconda $orthofinder/orthofinder.py -f $dir -t 16
