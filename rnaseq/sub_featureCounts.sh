#!/bin/bash
#Align RNAseq data with genome using featureCounts

#$ -S /bin/bash
#$ -cwd
#$ -pe smp 4
#$ -l virtual_free=1.2G
#$ -l h=blacklace02.blacklace|blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

Usage='qsub $ProgDir/sub_featureCounts.sh $BamFile $Gff $OutDir $OutFilePrefix'

# ---------------
# Step 1
# Collect inputs
# ---------------

InBam=$(basename $1)
InGff=$(basename $2)
OutDir=$3
Prefix=$4

CurDir=$PWD
WorkDir=$TMPDIR/featureCounts
echo "$WorkDir"
mkdir -p $WorkDir
cd $WorkDir

cp $1 ./
cp $2 ./


# ---------------
# Step 2
# Run featureCounts
# ---------------

featureCounts \
  -p -B -M -R \
  -T 4 \
  -a $InGff \
  -t exon \
  -g "Parent" \
  -o "$Prefix"_featurecounts.txt \
  $InBam

rm $InBam
rm $InGff
mkdir -p $CurDir/$OutDir
cp -r $WorkDir/* $CurDir/$OutDir/.
