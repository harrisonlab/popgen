#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=16G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace12.blacklace

##############################################
# Prep nanopore/pacbio mappings from bwa mem for SNP calling
### Remove multimapping reads
### add read group and sample name to each mapped read (preferably, the shortest ID possible)
#INPUT:
# 1st argument: input BAM file with your mappings
# 2nd argument: sample name (prefix) to be used to identify it in the future
# 3rd argument: output folder
#OUTPUT:
# Indexed BAM file with suffix "_rg" to be fed into SNP calling with GATK.
#############################################
input_bam=$1
prefix=$2

filename=$(basename "$input_bam")
name="${filename%.bam}"

### Output folder
cpath=$3
temp_dir="$TMPDIR"

### Prep
mkdir -p $temp_dir
cp $input_bam $temp_dir
cd $temp_dir

### Get rid of multimapping reads by filtering out on the XS:i: tag
samtools view -h $filename | awk '$17 !~ /XA:/|| $1 ~ /^@/' | samtools view -bS - > ${name}_unique.bam
samtools sort ${name}_unique.bam ${name}_unique_sorted
samtools index ${name}_unique_sorted.bam

### Add group and sample name (prefix)
picard=/home/sobczm/bin/picard-tools-2.5.0/picard.jar
java -jar $picard AddOrReplaceReadGroups \
INPUT=${name}_unique_sorted.bam \
OUTPUT=${name}_unique_sorted_rg.bam \
SORT_ORDER=coordinate CREATE_INDEX=true RGID=$prefix  RGSM=$prefix \
RGPL=Nanopore RGLB=library RGPU=barcode VALIDATION_STRINGENCY=LENIENT
samtools index ${name}_unique_sorted_rg.bam

### Cleanup
cp ${name}_unique_sorted_rg.bam $cpath
cp ${name}_unique_sorted_rg.bam.bai $cpath
rm -rf $temp_dir
