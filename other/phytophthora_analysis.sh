#!/bin/bash
scripts=/home/sobczm/bin/popgen/other
input=/home/sobczm/popgen/other/phytophthora

#Get OrthoMCL results into a more readable tabular format
python $scripts/rearrange_orthomcl_results_with_rubi.py All_Strains_plus_rubi_orthogroups.txt
python $scripts/rearrange_orthomcl_results_without_rubi.py All_Strains_orthogroups.txt
#How many orthogroups present in all strains but BC16?
Orthogroups=All_Strains_orthogroups.txt
cat $Orthogroups | grep 'A4|' | grep 'Bc1|' | grep 'Bc23|' | grep 'Nov27|' | grep 'Nov5|' | grep 'Nov71|' | grep 'Nov77|' | grep 'Nov9|' | grep 'ONT3|' | grep 'SCRP245_v2|' | grep -v 'Bc16|' >no_Bc16_groups.txt
#There appears to be 230 such orthogroups, out of 19952.
