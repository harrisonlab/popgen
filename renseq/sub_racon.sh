#!/bin/bash

# Alignment of minion reads to a minion assembly prior to running nanopolish variants
source /home/armita/.profile
PATH=${PATH}:/home/armita/prog/racon/racon/bin
PATH=${PATH}:/home/armita/prog/minimap/minimap
PATH=${PATH}:/home/armita/prog/minimap/minimap/utils/bin
PATH=${PATH}:/home/armita/prog/minimap/minimap/utils/plot
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l virtual_free=1G
#$ -l h=blacklace02.blacklace|blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace


Usage="sub_racon.sh <assembly.fa> <corrected_reads.fq.gz> <number_of_iterations> <output_directory>"
echo "$Usage"

# ---------------
# Step 1
# Collect inputs
# ---------------

AssemblyIn=$1
ReadsIn=$2
Iterations=$3
OutDir=$4

echo "Assembly - $AssemblyIn"
echo "Fasta reads - $ReadsIn"
echo "OutDir - $OutDir"

CurDir=$PWD
WorkDir=$TMPDIR/racon
mkdir -p $WorkDir
cd $WorkDir

Assembly=$(basename $AssemblyIn)
Assembly=$(echo $Assembly | sed 's/.utg/.fa/g')
Reads=$(basename $ReadsIn)
cp $CurDir/$AssemblyIn $Assembly
cp $CurDir/$ReadsIn $Reads

Prefix=$(echo $Assembly | cut -f1 -d '.')

mkdir -p $CurDir/$OutDir

cp $Assembly current-assembly.fa
for i in $(seq 1 $Iterations); do
  echo "Iteration - $i"
  minimap \
    -t16 \
    current-assembly.fa \
    $Reads \
    > racon_round_$i.reads_mapped.paf
  racon -t 16 $Reads racon_round_$i.reads_mapped.paf current-assembly.fa racon_round_$i.fasta
  cp racon_round_$i.fasta current-assembly.fa
  cp racon_round_$i.fasta $CurDir/$OutDir/"$Prefix"_racon_round_$i.fasta
done