#!/bin/bash
scripts=/home/sobczm/bin/popgen/other

cd /home/sobczm/popgen/other/rong
#Extract gene names to search for
for strain in 12008 51 53 58 61 JR2 Ls17 
do
cut -f1  /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/analysis/effectorP/V.dahliae/$strain/V.dahliae_${strain}_EffectorP_secreted_headers.txt | sed "s/^/$strain|/" >${strain}.txt
done

#Look for only orthogroups containing effector genes in at least one of the V. dahliae genomes
for a in 12008.txt 51.txt 53.txt 58.txt 61.txt JR2.txt Ls17.txt
do
grep -Fw -f $a All_Strains_orthogroups.txt >> effectors_orthogroups.txt
done

#Remove redundant identical entries
sort -k1 effectors_orthogroups.txt | uniq >effector_orthogroups_final.txt

#Change into Excel format
python $scripts/rearrange_orthomcl_results2.py effector_orthogroups_final.txt \
12008,,51,,53,,58,,61,,Ls17,,JR2,,VaMs102