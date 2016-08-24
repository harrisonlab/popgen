#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 16
#$ -l virtual_free=1G
#$ -l h=blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

#################################
# Stampy sensitive read mapping #
#################################
#			INPUT
assembly=$1
reads_for=$2
reads_rev=$3

##################################
#			OUTPUT				

In the current directory:

SAM file
sorted BAM file
BAI index
STAMPY log file
#################################


stampy=/home/sobczm/bin/stampy-1.0.29/stampy.py
assembly_file=$(basename "$assembly")
reads_for_file=$(basename "$reads_for")
reads_rev_file=$(basename "$reads_rev")

### Output folder
cpath=$PWD
temp_dir="$TMPDIR"

### Temp folder prep
mkdir -p $temp_dir
cp $assembly $reads_for $reads_rev $temp_dir
cd $temp_dir

### Build a genome file
$stampy -G "${assembly_file%.*}" $assembly_file

### Build a hash table
$stampy -g "${assembly_file%.*}" -H "${assembly_file%.*}"

### Map reads
$stampy -g "${assembly_file%.*}" -t 16 -h "${assembly_file%.*}" --substitutionrate=0.15 \
--xa-max=10 -f sam -o "${assembly_file%.*}_mapped_stampy.sam" -M $reads_for_file $reads_rev_file \
2>&1 | tee stampy_log.txt

### Get sorted BAM
samtools view -bS "${assembly_file%.*}_mapped_stampy.sam" -o "${assembly_file%.*}_mapped_stampy.bam"
samtools sort "${assembly_file%.*}_mapped_stampy.bam" "${assembly_file%.*}_mapped_stampy_sorted"
rm "${assembly_file%.*}_mapped_stampy.bam"
samtools index "${assembly_file%.*}_mapped_stampy_sorted.bam"

### Cleanup
mv *stampy* $cpath
rm -r $temp_dir
