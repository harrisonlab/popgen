#!/bin/bash

# Assemble PacBio data using SMRTdenovo

#$ -S /bin/bash
#$ -cwd
#$ -pe smp 16
#$ -l virtual_free=1.9G
#$ -l h=blacklace02.blacklace|blacklace06.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

source /home/armita/.profile
PATH=${PATH}:/home/armita/prog/smartdenovo/smartdenovo
Usage="sub_SMRTdenovo.sh <reads.fa> <outfile_prefix> <output_directory>"
echo "$Usage"

# ---------------
# Step 1
# Collect inputs
# ---------------

FastaIn=$1
Prefix=$2
OutDir=$3
echo  "Running SMARTdenovo with the following inputs:"
echo "FastaIn - $FastaIn"
echo "Prefix - $Prefix"
echo "OutDir - $OutDir"

CurPath=$PWD
WorkDir="$TMPDIR"/SMRTdenovo

# ---------------
# Step 2
# Run SMARTdenovo
# ---------------

mkdir -p $WorkDir
cd $WorkDir
Fasta=$(basename $FastaIn)
cp $CurPath/$FastaIn $Fasta

cat $Fasta | gunzip -cf > reads.fa
# fq2fa.pl reads.fq > reads.fa
smartdenovo.pl -t 16 reads.fa -p $Prefix > $Prefix.mak

make -f $Prefix.mak 2>&1 | tee "$Prefix"_run_log.txt

rm $Prefix.mak
rm $Fasta
rm reads.fa
mkdir -p $CurPath/$OutDir
for File in $(ls $WorkDir/wtasm*); do
  NewName=$(echo $File | sed "s/wtasm/$Prefix/g")
  mv $File $NewName
done
cp -r $WorkDir/* $CurPath/$OutDir/.