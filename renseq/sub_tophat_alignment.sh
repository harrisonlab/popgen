#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 4
#$ -l virtual_free=1G
#$ -l h=blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace


Usage="tophat_align.sh Genome.fa ReadFs.fq ReadRs.fq <Output_directory> [read_insert_gap_lgth] [stdev_read_insert_gap_lgth]"
Genome=$1
ReadF=$2
ReadR=$3
OutDir=$4

Organism=$(echo $ReadF | rev | cut -d "/" -f4 | rev)
Strain=$(echo $ReadF | rev | cut -d "/" -f3 | rev)

if [ $5 ] && [ $6 ]; then
  InsertGap=$5
  InsertStdDev=$6
  echo "$InsertGap"
  echo "$InsertStdDev"
  echo "running with an insert gap length (tophat -r option) and standard deviation of the\
  insert gap length (tophat --mate-std-dev toption) that was specified on the command \
  line."
else
  InsertGap='50'
  InsertStdDev='20'
  echo "$InsertGap"
  echo "$InsertStdDev"
  echo "running with tophat default insert gap length (tophat -r option) and \
  standard deviation of inset gap length (tophat --mate-std-dev toption). \
  It is suggested you run cufflinks on the aligned_bam file and use the \
  repoted Fragment Length Distribution and stdev to calculate these values \
  before resubmitting."
fi

OutName=$(basename "$Strain"_tophat)
echo "Internally storing in Dir= $OutName"

CurPath=$PWD
WorkDir=$TMPDIR/tophat_align

mkdir -p $WorkDir

echo "Copying genome file."
cp $Genome $WorkDir/"$OutName".fa

echo "Copying F reads"
cp $ReadF $WorkDir/ReadF

echo "Copying R reads"
cp $ReadR $WorkDir/ReadR

cd $WorkDir

echo "Unzip F read and R reads"
gunzip -cf ReadF > ReadF.fastq
gunzip -cf ReadR > ReadR.fastq

echo "Running bowtie"
bowtie2-build "$OutName".fa "$OutName"

echo "Running tophat"
tophat -o "$OutName" -p 4 --library-type fr-unstranded -r $InsertGap --mate-std-dev $InsertStdDev "$OutName" ReadF.fastq ReadR.fastq 2>&1 | tee $OutName/"$Strain"_report.txt
rm "$OutName".fa*
rm ReadF.fastq
rm ReadR.fastq

mkdir -p $OutDir
cp -r $WorkDir/$OutName/* $OutDir
rm -r $TMPDIR
