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
vvd=EAA28370
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
for gene in adv1 frq os4 vvd al1 sub1
do
cd $input/extended/$gene
python $scripts/prepare_gene_list.py $input/extended/$gene/${gene}.txt $gene
for a in *$gene; do sed -i '/^$/d' $a; done
#For each genome extract promoter sequences for given genes in the list.
for a in *$gene
do
fasta_file=$(echo $a | cut -d"." -f1)
python $scripts/keep_list_genes.py $a ../${fasta_file}*_promoters_2000.fasta
mv ../${fasta_file}*_promoters_2000_filtered.fasta ./
done
done

# Use FIMO to scan for motif GATCGA. Had to temporarily change p-value threshold,
# as the motif searched was very short, with the option --thresh 0.01
for gene in adv1 frq os4 vvd al1 sub1
do
cd $input/extended/$gene
for b in *_2000_filtered.fasta
do
qsub $scripts/sub_fimo.sh $input/extended/$gene/$b $gene GATCGA
done
done

#####Motif discovery - Top candidates and all candidates -> DREME 
#####Where are they found? (position)

#Establish the orthogroups containing those genes in each genome
dna=/home/sobczm/popgen/clock/DNA_genomes
pep=/home/sobczm/popgen/clock/pep_genomes
while read name;
do
grep "$name" $pep/OrthoFinder2/Results_Oct26/Orthogroups.csv
done <$dna/top_candidates.txt >top_candidates_orthogroups.txt

while read name;
do
grep "$name" $pep/OrthoFinder2/Results_Oct26/Orthogroups.csv
done <$dna/all_candidates.txt >all_candidates_orthogroups.txt

for gene in all_candidates top_candidates
do
cd $input/extended/$gene
python $scripts/prepare_gene_list.py $input/extended/$gene/${gene}.txt $gene
for a in *$gene; do sed -i '/^$/d' $a; done
#For each genome extract promoter sequences for given genes in the list.
for a in *.$gene
do
echo $a
fasta_file=$(echo $a | cut -d"." -f1)
echo $fasta_file
python $scripts/keep_list_genes.py $a ../${fasta_file}*_promoters_2000.fasta
mv ../${fasta_file}*_promoters_2000_filtered.fasta ./
done
done

#Count the number of sequences in each file and create an appropriate background
#control file with the same number of random promoters.
for file in $input/extended/all_candidates/*filtered.fasta 
do
lines=`wc -l $file | cut -f1 -d' '`
ngenes="$((lines/2))"
ffilename=$(basename "$file")
assembly=${ffilename%_filtered.fasta}.fasta
qsub $scripts/sub_fasta_subsample.sh $input/extended/$assembly $ngenes
done

for file in $input/extended/top_candidates/*filtered.fasta 
do
lines=`wc -l $file | cut -f1 -d' '`
ngenes="$((lines/2))"
ffilename=$(basename "$file")
assembly=${ffilename%_filtered.fasta}.fasta
qsub $scripts/sub_fasta_subsample.sh $input/extended/$assembly $ngenes
done

#DREME motif analysis
for control in $input/extended/all_candidates/*random*.fasta
do
file=${control%_random*}_filtered.fasta
qsub $scripts/sub_dreme.sh $file $control
done

for control in $input/extended/top_candidates/*random*.fasta
do
file=${control%_random*}_filtered.fasta
qsub $scripts/sub_dreme.sh $file $control
done
