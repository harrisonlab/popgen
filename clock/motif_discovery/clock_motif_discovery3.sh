#!/bin/bash
input=/home/sobczm/popgen/clock/DNA_genomes/promoters
scripts=/home/sobczm/bin/popgen/clock/motif_discovery
meme=/home/sobczm/bin/meme_4.11.2/bin

#Look for motif GATC in the promoters from the first part of the file
#clock_motif_discovery2.sh plus the following genes:
#eas - NCU08457
#cry - NCU00582
#con-8 - NCU09235
#bli-4 - NCU08699

#Extract the promoters of each gene in each species
cd $input/extended
#Find the protein IDS of each gene. 
eas=EAA34064
cry=EAA36486
con8=EAA34576
bli4=EAA34324

#Establish the orthogroups containing those genes in each genome
pep=/home/sobczm/popgen/clock/pep_genomes
head -1 $pep/OrthoFinder2/Results_Oct26/Orthogroups.csv >>eas.txt
grep $eas $pep/OrthoFinder2/Results_Oct26/Orthogroups.csv >>eas.txt
head -1 $pep/OrthoFinder2/Results_Oct26/Orthogroups.csv >>cry.txt
grep $cry $pep/OrthoFinder2/Results_Oct26/Orthogroups.csv >>cry.txt
head -1 $pep/OrthoFinder2/Results_Oct26/Orthogroups.csv >>con8.txt
grep $con8 $pep/OrthoFinder2/Results_Oct26/Orthogroups.csv >>con8.txt
head -1 $pep/OrthoFinder2/Results_Oct26/Orthogroups.csv >>bli4.txt
grep $bli4 $pep/OrthoFinder2/Results_Oct26/Orthogroups.csv >>bli4.txt

#Fish out the homologs of the genes in other species
for gene in eas cry con8 bli4
do
cd $input/extended/$gene
python $scripts/prepare_gene_list.py $input/extended/$gene/${gene} $gene
for a in *$gene; do sed -i '/^$/d' $a; done
#For each genome extract promoter sequences for given genes in the list.
for a in *$gene
do
fasta_file=$(echo $a | cut -d"." -f1)
python $scripts/keep_list_genes.py $a ../${fasta_file}*_promoters_2000.fasta
mv ../${fasta_file}*_promoters_2000_filtered.fasta ./
done
done

# Use FIMO to scan for motif GATCGA. Had to temporarily change p-value threshold, as the motif searched was very short, with the option --thresh 0.1
#Also, searched only on plus strand as the sequence is palindromic,
for gene in adv1 frq os4 vvd al1 sub1 eas cry con8 bli4
do
cd $input/extended/$gene
for b in *_2000_filtered.fasta
do
qsub $scripts/sub_fimo.sh $input/extended/$gene/$b LRE GATC
done
done

#Use Glam2Scan to look for extended gap Cbox motif:
#CGAT-----CCGCT
#spacing between 0-30 bp
cd $input/extended/frq
for a in *_promoters_2000_filtered.fasta; do qsub $scripts/sub_glam2scan.sh $a ../cbox_extended.txt; done
cd $input/extended/vvd
for a in *_promoters_2000_filtered.fasta; do qsub $scripts/sub_glam2scan.sh $a ../cbox_extended.txt; done

#Use Glam2Scan to look for extended gap motif from He (2005):
#GATNC-----CGATN
#spacing between 0-30 bp
cd $input/extended/frq
for a in *_promoters_2000_filtered.fasta; do qsub $scripts/sub_glam2scan.sh $a ../he2005_extended.txt; done
cd $input/extended/vvd
for a in *_promoters_2000_filtered.fasta; do qsub $scripts/sub_glam2scan.sh $a ../he2005_extended.txt; done