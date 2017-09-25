#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 3
#$ -l virtual_free=20G
#$ -l h=blacklace02.blacklace|blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

lorma=/home/sobczm/bin/LoRMA-0.4
input_fasta=$1
OutDir=$2

CurPath=$PWD
mkdir -p $WorkDir
cd $WorkDir
Fasta=$(basename $input_fasta)
cp $CurPath/$input_fasta $Fasta

$lorma/lorma.sh -threads 3 $Fasta

mkdir -p $CurPath/$OutDir
cp $WorkDir/* $CurPath/$OutDir