#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l virtual_free=1.9G
#$ -l h=blacklace02.blacklace|blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace11.blacklace

source /home/armita/.profile
Usage="sub_canu_renseq.sh <reads.fq> <Genome_size[e.g.45m]> <outfile_prefix> <output_directory> [<specification_file.txt>]"
echo "$Usage"

# ---------------
# Step 1
# Collect inputs
# ---------------
canu=/home/armita/prog/canu/canu-1.6/Linux-amd64/bin/canu
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

cpath=$PWD
mkdir $OutDir
temp_dir="$TMPDIR"
mkdir -p $temp_dir

cp $FastqIn $temp_dir

Fastq=$(basename "$FastqIn")
cd $temp_dir

# ---------------
# Step 2
# Run Canu
# ---------------

$canu \
  -correct \
  -useGrid=false \
  -overlapper=mhap \
  -utgReAlign=true \
  -d $OutDir \
  -p $Prefix \
  genomeSize="$Size" \
  -nanopore-raw $Fastq \
  2>&1 | tee canu_run_log.txt

  $canu \
  -trim \
  -useGrid=false \
  -overlapper=mhap \
  -utgReAlign=true \
  -d $OutDir \
  -p $Prefix \
  genomeSize="$Size" \
  -nanopore-corrected $OutDir/$Prefix.correctedReads.fasta.gz \
  2>&1 | tee canu_run_log.txt

  rm *.fastq
  cp -r $OutDir/* $cpath/$OutDir
  rm -rf $temp_dir
