#!/bin/bash
input=/home/sobczm/popgen/clock/DNA_genomes/promoters
scripts=/home/sobczm/bin/popgen/clock/motif_discovery
meme=/home/sobczm/bin/meme_4.11.2/bin

#Search within 2 kbp of start site.
############Motif scanning - a couple of key ones like:
#adv1 (NCU07392), os-4 (NCU03071), frq (NCU02265), vvd (NCU03967), 
#al-1 (NCU00552), sub-1 (NCU01154)
#GATCGA => its position varied.
#Extract the promoters of each gene in each species
cd $input/extended
#Find the protein IDS of each gene. 
adv1=EAA29240
os4=EAA34928
frq=ESA42015
vvd=ESA42013
al1=EAA35477
sub1=ESA42508

#Establish the orthogroups containing those genes in each genome
pep=/home/sobczm/popgen/clock/pep_genomes
head -1 $pep/OrthoFinder2/Results_Oct26/Orthogroups.csv >>adv1.txt
grep $adv1 $pep/OrthoFinder2/Results_Oct26/Orthogroups.csv >>adv1.txt
head -1 $pep/OrthoFinder2/Results_Oct26/Orthogroups.csv >>os4.txt
grep $os4 $pep/OrthoFinder2/Results_Oct26/Orthogroups.csv >>os4.txt
head -1 $pep/OrthoFinder2/Results_Oct26/Orthogroups.csv >>frq.txt
grep $frq $pep/OrthoFinder2/Results_Oct26/Orthogroups.csv >>frq.txt
head -1 $pep/OrthoFinder2/Results_Oct26/Orthogroups.csv >>vvd.txt
grep $vvd $pep/OrthoFinder2/Results_Oct26/Orthogroups.csv >>vvd.txt
head -1 $pep/OrthoFinder2/Results_Oct26/Orthogroups.csv >>al1.txt
grep $al1 $pep/OrthoFinder2/Results_Oct26/Orthogroups.csv >>al1.txt
head -1 $pep/OrthoFinder2/Results_Oct26/Orthogroups.csv >>sub1.txt
grep $sub1 $pep/OrthoFinder2/Results_Oct26/Orthogroups.csv >>sub1.txt


#Fish out the homologs of the genes in other species
#To accomplish that, first parse the input spreadsheet with protein names
#to produce a custom list for each species.
python $scripts/prepare_gene_list.py $input/extended/frq/frq.txt frq
for a in *frq; do sed -i '/^$/d' $a; done
#For each genome extract promoter sequences for given genes in the list.
for a in *frq
do
fasta_file=$(echo )
python $scripts/keep_list_genes.py $input/ace/Botrytis_cinerea.pep.fa.ace \
Botrytis_cinerea.ASM83294v1.dna_rm.toplevel_promoters_2000.fasta
done



#If any new variations on the motif discovered here, check for
#motif enrichment in top candidates and all candidates.



#Motif discovery - Top candidates and all candidates -> DREME 
#Where are they found? (position)


