#!/bin/bash
input=/home/sobczm/popgen/clock/coding_genomes
scripts=/home/sobczm/bin/popgen/clock

####Calculation of pairwise Ka/Ks ratios between V. dahliae and V. alfalfae (albo-atrum)
####genes involved in microslerotia development

#Find the orthology groups containing the relevant genes
pep=/home/sobczm/popgen/clock/pep_genomes
while read name;
do
grep "$name" $pep/OrthoFinder2/Results_Oct26/Orthogroups.csv
done <$input/microsclerotia.txt >microsclerotia_orthogroups.txt

# Verticillum alfalfae and dahliae are columns 28 and 29 respectively,
# and rows are individual orthogroups.

### Prepare input FASTA for pairwise Ka/Ks calculations. 
#Independent FASTA file for all possible pairwise (ie. one sequence per species) comparisons
#in an orthogroup.

python $scripts/dn_ds/pairwise_ka_ks.py --o $input/microsclerotia_orthogroups.txt --col1 28 --fasta1 $input/Verticillium_alfalfae_cds.fasta --prefix1 Va --col2 29 --fasta2 $input/Verticillium_dahliae_cds.fasta --prefix2 Vd

