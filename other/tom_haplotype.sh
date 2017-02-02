#!/bin/bash

input_vcf=/home/sobczm/popgen/other/phytophthora/95m_contigs_unmasked_UK123_filtered.recode_annotated_haplo.vcf
input_genome=/home/sobczm/popgen/other/phytophthora/genomes/95m_contigs_hardmasked.fa
scripts=/home/sobczm/bin/popgen/summary_stats

cd /home/sobczm/popgen/other/phytophthora
mkdir 4GT && cd 4GT

#Create PopGenome input 
python $scripts/vcf_to_fasta.py $input_vcf $input_genome 2

mkdir contigs && mv *.fasta ./contigs
cd contigs
for f in *.fasta
do
folder=${f%.fasta}
mkdir $folder
mv $f $folder
done
#Gff files
cd ..
gff=/home/groups/harrisonlab/project_files/phytophthora_fragariae/gene_pred/codingquary/P.fragariae/Bc16/final/final_genes_appended.gff3
$scripts/split_gff_contig.sh $gff
mkdir gff && mv *.gff ./gff
#Check for orphan contigs with no matching gff file, which need to be removed prior to the run.
for a in $PWD/contigs/*/*.fasta
do
filename=$(basename "$a")
expected_gff="$PWD/gff/${filename%.fa*}.gff"
if [ ! -f "$expected_gff" ];
then
   rm -rf $(dirname $a)
fi
done

#Run a script to calculate the four gamete test
