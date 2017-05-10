#!/bin/bash
input=/home/sobczm/popgen/rnaseq
scripts=/home/sobczm/bin/popgen/rnaseq
cd $input
#Count up the reads using featureCounts
Gff=$input/Fragaria_vesca_v1.1.a2.gff3
for BamFile in $(ls $input/pcac_*/vesca/star_aligmentAligned.sortedByCoord.out.bam); do
  OutDir=$(dirname $BamFile)
  OutDir=$OutDir/featureCounts
  Prefix=$(echo $BamFile | rev | cut -f3 -d '/' | rev)
  echo $OutDir
  echo $Prefix
  qsub $scripts/sub_featureCounts.sh $BamFile $Gff $OutDir $Prefix
done

#A file was created with columns referring to experimental treatments:
OutDir=alignment/star/P.cactorum/414_v2/DeSeq
mkdir -p $OutDir
# make file in excel and copy accross
cat $OutDir/P.cactorum_RNAseq_design.txt | tr -d '\r' | sed 's/Timepoint/Timepoint\n/g' | sed "s/hours/hours\n/g" > $OutDir/P.cactorum_RNAseq_design_parsed.txt
for File in $(ls alignment/star/P.cactorum/414_v2/*/*/featureCounts/*_featurecounts.txt); do
  echo $File;
  cp $CurDir/$File $OutDir/.;
done
for File in $(ls $OutDir/*_featurecounts.txt); do
  Prefix=$(echo $File | rev | cut -f1 -d '/' | rev | sed 's/_featurecounts.txt//g' | sed "s/_totRNA_S.*_L/_L/g")
  sed -ie "s/star_aligmentAligned.sortedByCoord.out.bam/$Prefix/g" $File
done