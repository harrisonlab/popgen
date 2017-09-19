#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 16
#$ -l virtual_free=1G
#$ -l h=blacklace02.blacklace|blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace


# Wrapper to submit Porechop jobs to the cluster
# Demultiplex and trim adapters
# Porechop splits MinION reads at locations of adapter sequences

Usage="sub_porechop.sh <raw_reads.fastq.gz> <output_directory>"
echo $Usage

#######  Step 1	 ########
# Initialise values
#########################
porechop=/home/armita/prog/porechop/Porechop
RawReads=$(basename $1)
OutDir=$2
CurDir=$PWD
WorkDir=$TMPDIR/porechop

Prefix=$(echo $RawReads | cut -f1 -d '.')

echo "Raw reads - $RawReads"
echo "Output directory - $OutDir"
echo "Output files will carry the prefix - $Prefix"

#######  Step 2	 ########
# Copy data on worker node
#########################

mkdir -p $WorkDir
cp $1 $WorkDir
cd $WorkDir
#######  Step 3	 ########
# Run Porechop
#########################

$porechop/porechop-runner.py -i $RawReads -b $OutDir --threads 16 --adapter_threshold 80.0 > "$Prefix"_trim_log.txt

#######  Step 3	 ########
# Tidy up
#########################
cp -r $OutDir $CurDir
cp "$Prefix"_trim_log.txt $CurDir/$OutDir