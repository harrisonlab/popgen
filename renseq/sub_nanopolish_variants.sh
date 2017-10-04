#!/bin/bash

# Alignment of minion reads to a minion assembly prior to running nanopolish variants

#$ -S /bin/bash
#$ -cwd
#$ -pe smp 4
#$ -l virtual_free=1G
# wont submit to node 1, 2, 11
#$ -l h=blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace
# wont submit to node 1, 11
# #$ -l h=blacklace02.blacklace|blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

Usage="sub_nanopolish_variants.sh <assembly.fa> <reads.fa.gz> <aligned_ONTreads.sam> <ploidy[e.g.1]> <region> <output_directory>"
echo "$Usage"

# ---------------
# Step 1
# Collect inputs
# ---------------

AssemblyIn=$1
ReadsIn=$2
AlignedIn=$3
Ploidy=$4
Region=$5
OutDir=$6
Prefix=$(echo $Region | sed -e "s/ /_/g")

echo "Assembly - $AssemblyIn"
echo "Fasta reads - $ReadsIn"
echo "Aligned reads - $AlignedIn"
echo "OutDir - $OutDir"



CurDir=$PWD
WorkDir=$TMPDIR/nanopolish
# WorkDir=/tmp/nanopolish
mkdir -p $WorkDir
cd $WorkDir

Assembly=$(basename $AssemblyIn)
Reads=$(basename $ReadsIn)
Aligned=$(basename $AlignedIn)
cp $CurDir/$AssemblyIn $Assembly
cp $CurDir/$AlignedIn $Aligned
cp $CurDir/$AlignedIn.bai $Aligned.bai
cp $CurDir/$ReadsIn $Reads

cat $Reads | gunzip -cf > reads.fa
source /home/sobczm/.generic_profile
which samtools
which nanopolish
#nanopolish=/home/armita/prog/nanopolish/nanopolish_edited/nanopolish
nanopolish index reads.fa
nanopolish variants \
  -t 4 \
  --reads reads.fa \
  --snps \
  --bam $Aligned \
  --genome $Assembly \
  --ploidy $Ploidy \
  -w $Region \
  --fix-homopolymers \
  --min-candidate-frequency 0.1 \
  --min-candidate-depth 10 \
  --max-haplotypes 1000000 \
  --calculate-all-support > "$Prefix"_variants.txt

mkdir -p $CurDir/$OutDir
cp "$Prefix"* $CurDir/$OutDir/.
