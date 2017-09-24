#!/bin/bash

# Alignment of minion reads to a minion assembly prior to running nanopolish variants

#$ -S /bin/bash
#$ -cwd
#$ -pe smp 16
#$ -l virtual_free=1G
#$ -l h=blacklace02.blacklace|blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

source /home/armita/.profile

Usage="sub_bwa_nanopolish.sh <assembly.fa> <reads.fa.gz> <output_directory>"
echo "$Usage"

# ---------------
# Step 1
# Collect inputs
# ---------------

AssemblyIn=$1
ReadsIn=$2
OutDir=$3

echo "Assembly - $AssemblyIn"
echo "Fasta reads - $ReadsIn"
echo "OutDir - $OutDir"

CurDir=$PWD
WorkDir=$TMPDIR/nanopolish
mkdir -p $WorkDir
cd $WorkDir

Assembly=$(basename $AssemblyIn)
Reads=$(basename $ReadsIn)
cp $CurDir/$AssemblyIn $Assembly
cp $CurDir/$ReadsIn $Reads

mkdir -p $CurDir/$OutDir

bwa index $Assembly
bwa mem -x ont2d -t 16 $Assembly $Reads > aligned.sam
samtools view -b -S aligned.sam | samtools sort - -o tmp > reads.sorted.bam
samtools index reads.sorted.bam

rm $Assembly
rm $Reads
rm aligned.sam

cp -r $WorkDir/* $CurDir/$OutDir/.
