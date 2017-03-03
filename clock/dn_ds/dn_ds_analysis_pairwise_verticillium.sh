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

#Move all the output FASTA file into a new folder
mkdir $input/microsclerotia
mv O*.fasta $input/microsclerotia

#Remove terminal stop codons
for a in *.fasta
do
python $scripts/dn_ds/remove_terminal_stop.py $a
done

#Transalign the sequenes by codon
for a in *ns.fasta
do
sh $scripts/dn_ds/transalign.sh $a 
done
#has to be extended phylip format for codeml

#Write codeml control files (calculate omega and omega fixed)
perl $scripts/dn_ds/write_codeml_control.pl

#Run codeml
codeml=/home/sobczm/bin/paml4.8/bin/codeml
for a in *.ctl
do
$codeml $PWD/$a
done

#Parse codeml output
perl $scripts/dn_ds/parse_codeml_result.pl
