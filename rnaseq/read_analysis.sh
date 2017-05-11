#!/bin/bash
input=/home/sobczm/popgen/rnaseq
scripts=/home/sobczm/bin/popgen/rnaseq
cd $input
#Count up the reads using featureCounts
Gff=$input/Fragaria_vesca_v1.1.a2.gff3
BamFile=star_aligmentAligned.sortedByCoord.out.bam
for k in ./pcac_*/vesca
do
cd $PWD/$k
Prefix=$(echo $PWD/$BamFile | rev | cut -f3 -d '/' | rev)
echo $Prefix
qsub $scripts/sub_featureCounts.sh $PWD/$BamFile $Gff $Prefix
cd $input
done

Alignment stats were collected:
#First round of alignment
cd $input
for File in $(ls vesca_*/star_aligmentLog.final.out); do
Sample=$(echo $File | rev | cut -f2 -d '/' | rev | sed 's/vesca_//');
InputReads=$(cat $File | grep 'Number of input reads' | cut -f2);
ReadNumU=$(cat $File | grep 'Uniquely' | grep 'number' | cut -f2);
ReadPercU=$(cat $File | grep 'Uniquely' | grep '%' | cut -f2);
ReadNumM=$(cat $File | grep 'multiple' | grep 'Number' | cut -f2);
ReadPercM=$(cat $File | grep 'multiple' | grep '%' | cut -f2);
Mismatch=$(cat $File | grep 'Mismatch rate per base' | grep '%' | cut -f2);
echo -e "$Sample""\t""$InputReads""\t" "$ReadNumU""\t""$ReadPercU""\t""$ReadNumM""\t""$ReadPercM""\t""$Mismatch";  
done
#Second round of alignment
cd $input
for File in $(ls pcac_*/vesca/star_aligmentLog.final.out); do
Sample=$(echo $File | rev | cut -f3 -d '/' | rev | sed 's/pcac_//');
InputReads=$(cat $File | grep 'Number of input reads' | cut -f2);
ReadNumU=$(cat $File | grep 'Uniquely' | grep 'number' | cut -f2);
ReadPercU=$(cat $File | grep 'Uniquely' | grep '%' | cut -f2);
ReadNumM=$(cat $File | grep 'multiple' | grep 'Number' | cut -f2);
ReadPercM=$(cat $File | grep 'multiple' | grep '%' | cut -f2);
Mismatch=$(cat $File | grep 'Mismatch rate per base' | grep '%' | cut -f2);
echo -e "$Sample""\t""$InputReads""\t" "$ReadNumU""\t""$ReadPercU""\t""$ReadNumM""\t""$ReadPercM""\t""$Mismatch";  
done

#Copy featureCounts tables into one directory (second round of alignment) and rename to remove pcac_ prefix
for k in ./pcac_*/vesca/*_featurecounts.txt ./pcac_*/vesca/*_featurecounts.txt.summary
do
cp $k ./featureCounts
fb=$(basename $k)
new_name=$(echo $fb | sed 's/pcac_//')
mv ./featureCounts/$fb ./featureCounts/$new_name
done

#A file was created with columns referring to experimental treatments. Copy it over
cd $input/featureCounts
OutDir=/home/groups/harrisonlab/project_files/idris/alignment/star/P.cactorum/414_v2/DeSeq
cp $OutDir/P.cactorum_RNAseq_design_parsed.txt $input/featureCounts

for a in *featurecounts.txt
do
cp $a ${a}_original
done

for File in $(ls *_featurecounts.txt)
do
  Prefix=$(echo $File | rev | cut -f1 -d '/' | rev | sed 's/_featurecounts.txt//g' | sed "s/_totRNA_S.*_L/_L/g")
  echo $Prefix
  sed -ie "s/star_aligmentAligned.sortedByCoord.out.bam/$Prefix/g" $File
  sed -i "s/pcac_//g" $File
done


##Carry out the intial DeSeq analysis with script DeSeq_vesca.R