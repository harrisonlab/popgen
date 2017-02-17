#$ -S /bin/bash
#$ -cwd
#$ -pe smp 16
#$ -l virtual_free=1.9G
#$ -l h=blacklace11.blacklace

Usage="sub_canu_renseq.sh <reads.fq> <Genome_size[e.g.45m]> <outfile_prefix> <output_directory> [<specification_file.txt>]"
echo "$Usage"

# ---------------
# Step 1
# Collect inputs
# ---------------
canu=/home/armita/prog/canu/canu-1.3/Linux-amd64/bin/canu

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


# ---------------
# Step 2
# Run Canu
# ---------------

$canu \
  useGrid=false \
  $AdditionalCommands \
  -d $OutDir \
  -p $Prefix genomeSize="$Size" \
  -pacbio-raw $FastqIn \
  2>&1 | tee canu_run_log.txt
