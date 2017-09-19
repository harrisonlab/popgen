#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l virtual_free=10G
#$ -l h=blacklace02.blacklace|blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

##Nanocorrect for slow self-correction of Nanopore reads.

Usage="sub_nanocorrect.sh <raw_reads.fasta> <output_directory_prefix>"
echo $Usage

nanocorrect=/home/sobczm/bin/nanocorrect
RawReads=$(basename $1)
Prefix=$2
CurDir=$PWD
WorkDir=$TMPDIR/nanocorrect

mkdir -p $WorkDir
cp $1 $WorkDir
cd $WorkDir

$nanocorrect/nanocorrect-overlap INPUT=$RawReads NAME=$Prefix
$nanocorrect/nanocorrect.py nc $Prefix > ${RawReads%.fasta}_nanocorrect.fasta

cp ${RawReads%.fasta}_nanocorrect.fasta $CurDir
