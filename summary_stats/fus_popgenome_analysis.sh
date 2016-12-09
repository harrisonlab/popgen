#!/bin/bash
input=/home/sobczm/popgen/summary_stats/noA13
scripts=/home/sobczm/bin/popgen/summary_stats

# In order to calculate different statistics in Popgenome, the input has to be arranged
# in a particular way.
# The input directory should contain two folders.
# Folder No. 1: named "gff", contains GFF files for all the contigs output from the split_gff_contig.sh script
# Folder No. 2: named "contigs", contains subfolders, each subfolder named with exact contig name
# and containing one individual contig FASTA file, also named with exact contig name, as output from vcf_to_fasta.py

#An example on how to create this directory structure

cd $input/all
#This folder contaings only contig FASTA files
#So create a new "contigs" directory to hold those files:
mkdir contigs
mv *.fasta ./contigs

#copy the "gff" folder containing gff files
cp -r /home/sobczm/popgen/summary_stats/gff ./

#The last step: in the folder "contigs" create subfolders, each to hold one contig FASTA file
cd contigs
for f in *.fasta
do
folder=${f%.fasta}
mkdir $folder
mv $f $folder
done

#Navigate to the input folder holding the two folders: "contigs" and "gff"
#to proceed with Popgenome run.
cd $input/all
#The R script used below is custom-made for each run (see first few lines of it)
#It requires custom definition of populations, and individual assignment to them.
#The example below calculates nucleotide diversity within (Pi) and between (Dxy) populations/
#Other scripts (sub_calculate_neutrality_stats.sh) are used in analogous manner.

qsub $scripts/sub_calculate_nucleotide_diversity.sh
qsub $scripts/sub_calculate_neutrality_stats.sh
qsub $scripts/sub_calculate_fst.sh
qsub $scripts/sub_calculate_haplotype_based_stats.sh

#This calculation was done over all sites. Now going to proceed for site subsets:
#synonymous, non-synonymous and four-fold degenerate (silent), in the respective folders

#four_fold_degenerate (analogous to above, for all sites)
cd $input/silent
mkdir contigs
mv *.fasta ./contigs
cp -r /home/sobczm/popgen/summary_stats/gff ./
cd contigs
for f in *.fasta
do
folder=${f%.fasta}
mkdir $folder
mv $f $folder
done
cd $input/silent

qsub $scripts/sub_calculate_nucleotide_diversity.sh
qsub $scripts/sub_calculate_neutrality_stats.sh
qsub $scripts/sub_calculate_fst.sh
qsub $scripts/sub_calculate_haplotype_based_stats.sh

#For synonymous and non-synonymous have to create FASTA input first, as done
#for silent and all sites in fus_variant_annotation.sh
##synonymous
cd $input
ref_genome=/home/sobczm/popgen/summary_stats/Fus2_canu_contigs_unmasked.fa
python $scripts/vcf_to_fasta.py Fus2_canu_contigs_unmasked_noA13_filtered.recode_annotated_syn.vcf $ref_genome 1
#Moving each subset of FASTA files into a separate dir.
mkdir syn
mv *.fasta ./syn

##non-synonymous
cd $input
ref_genome=/home/sobczm/popgen/summary_stats/Fus2_canu_contigs_unmasked.fa
python $scripts/vcf_to_fasta.py Fus2_canu_contigs_unmasked_noA13_filtered.recode_annotated_nonsyn.vcf $ref_genome 1
#Moving each subset of FASTA files into a separate dir.
mkdir nonsyn
mv *.fasta ./nonsyn

## And now back to creating dir structure and carrying out Popgenome analysis
cd $input/syn
mkdir contigs
mv *.fasta ./contigs
cp -r /home/sobczm/popgen/summary_stats/gff ./
cd contigs
for f in *.fasta
do
folder=${f%.fasta}
mkdir $folder
mv $f $folder
done

cd $input/syn
qsub $scripts/sub_calculate_nucleotide_diversity.sh
qsub $scripts/sub_calculate_neutrality_stats.sh
qsub $scripts/sub_calculate_fst.sh
qsub $scripts/sub_calculate_haplotype_based_stats.sh

cd $input/nonsyn
mkdir contigs
mv *.fasta ./contigs
cp -r /home/sobczm/popgen/summary_stats/gff ./
cd contigs
for f in *.fasta
do
folder=${f%.fasta}
mkdir $folder
mv $f $folder
done

cd $input/nonsyn
qsub $scripts/sub_calculate_nucleotide_diversity.sh
qsub $scripts/sub_calculate_neutrality_stats.sh
qsub $scripts/sub_calculate_fst.sh
qsub $scripts/sub_calculate_haplotype_based_stats.sh
