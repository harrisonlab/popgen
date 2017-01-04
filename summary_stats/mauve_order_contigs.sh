#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l virtual_free=1G
#$ -l h=blacklace02.blacklace|blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

MauveDir=$1
Reference=$2
Genome=$3
OutDir=$4

echo $MauveDir
echo $Reference
echo $Genome
echo $OutDir

LogName=$(basename $Genome)

cd $MauveDir
mkdir -p $OutDir
java -Xmx500m -Djava.awt.headless=true -cp Mauve.jar org.gel.mauve.contigs.ContigOrderer -output $OutDir -ref $Reference -draft $Genome  2>&1 $OutDir/$LogName.order_config.log
