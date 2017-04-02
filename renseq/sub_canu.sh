#!/bin/bash

# Assemble PacBio data using Canu

#$ -S /bin/bash
#$ -cwd
#$ -pe smp 16
#$ -l virtual_free=1.9G
#$ -l h=blacklace02.blacklace|blacklace06.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace


Usage="submit_canu.sh <reads.fq> <Genome_size[e.g.45m]> <outfile_prefix> <output_directory> [<specification_file.txt>]"
echo "$Usage"

# ---------------
# Step 1
# Collect inputs
# ---------------

FastqIn=$1
Size=$2
Prefix=$3
OutDir=$4
AdditionalCommands=""
if [ $5 ]; then
  SpecFile=$5
  AdditionalCommands="-s $SpecFile"
fi
echo  "Running Canu with the following inputs:"
echo "FastqIn - $FastqIn"
echo "Size - $Size"
echo "Prefix - $Prefix"
echo "OutDir - $OutDir"

CurPath=$PWD
WorkDir="$TMPDIR"/canu

# ---------------
# Step 2
# Run Canu
# ---------------

mkdir -p $WorkDir
cd $WorkDir
Fastq=$(basename $FastqIn)
cp $CurPath/$FastqIn $WorkDir/$Fastq

canu \
  useGrid=false \
  $AdditionalCommands \
  -d $WorkDir/assembly \
  -p $Prefix genomeSize="$Size" \
  -pacbio-raw $Fastq \
  2>&1 | tee canu_run_log.txt

mkdir -p $CurPath/$OutDir
cp canu_run_log.txt $CurPath/$OutDir/.
cp $WorkDir/assembly/* $CurPath/$OutDir/.