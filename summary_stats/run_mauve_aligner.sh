#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l virtual_free=30G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace12.blacklace

mauve=/home/armita/prog/mauve/mauve_snapshot_2015-02-13/linux-x64

OutDir=$1
GenomeList=$2

mkdir -p $OutDir

ScratchPath1=$TMPDIR/tmp1
ScratchPath2=$TMPDIR/tmp2
mkdir -p $ScratchPath1
mkdir -p $ScratchPath2

$mauve/mauveAligner --weight=300 --output=$OutDir/output.mauve --output-alignment=$OutDir/aligned_genomes.xmfa \
--scratch-path $ScratchPath1 --scratch-path $ScratchPath2 $GenomeList
