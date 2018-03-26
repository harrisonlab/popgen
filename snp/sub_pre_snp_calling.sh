#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=16G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace12.blacklace

##############################################
# Prep mappings from Bowtie2 for SNP calling
### Remove multimapping reads, discordant reads. PCR and optical duplicates, and
### add read group and sample name to each mapped read (preferably, the shortest ID possible)
#INPUT:
# 1st argument: input SAM file with your mappings
# 2nd argument: sample name (prefix) to be used to identify it in the future
#OUTPUT:
# Indexed BAM file with suffix "nodup_rg" to be fed into SNP calling with GATK.
#############################################
source /home/sobczm/.profile
input_sam=$1
prefix=$2
filename=$(basename "$input_sam")
name="${filename%.*}"

### Output folder
cpath=$PWD
temp_dir="$TMPDIR"

### Prep
mkdir -p $temp_dir
cp $input_sam $temp_dir
cd $temp_dir

### If the provided file is a .bam file then convert it into sam format.

extension=$(echo $filename | rev | cut -f1 -d '.' | rev)
if [[ $Extension == "bam" ]]; then
  samtools view -o $name.sam $filename
  filename=$name.sam
fi

### Get rid of multimapping reads by filtering out on the XS:i: tag
grep -v "XS:i" $filename >temp && mv temp $filename
samtools view -bS -o $name.bam $filename
samtools sort $name.bam ${name}_nomulti_sorted
samtools index ${name}_nomulti_sorted.bam

### Keep only reads with "paired reads" and "properly paired reads" flags.
samtools view -b -h -f 3 -o ${name}_nomulti_proper.bam ${name}_nomulti_sorted.bam
### Sort for downstream analyses
samtools sort ${name}_nomulti_proper.bam ${name}_nomulti_proper_sorted
samtools index ${name}_nomulti_proper_sorted.bam

### Remove PCR and optical duplicates
picard=/home/sobczm/bin/picard-tools-2.5.0/picard.jar
java -jar $picard MarkDuplicates \
INPUT=${name}_nomulti_proper_sorted.bam
OUTPUT=${name}_nomulti_proper_sorted_nodup.bam \
METRICS_FILE=${name}_nomulti_proper_sorted_nodup.txt \
REMOVE_DUPLICATES=TRUE ASSUME_SORTED=TRUE MAX_RECORDS_IN_RAM=500000000 \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 VALIDATION_STRINGENCY=LENIENT
### Add group and sample name (prefix)
java -jar $picard AddOrReplaceReadGroups \
INPUT=${name}_nomulti_proper_sorted_nodup.bam \
OUTPUT=${name}_nomulti_proper_sorted_nodup_rg.bam \
SORT_ORDER=coordinate CREATE_INDEX=true RGID=$prefix  RGSM=$prefix \
RGPL=Illumina RGLB=library RGPU=barcode VALIDATION_STRINGENCY=LENIENT
samtools index ${name}_nomulti_proper_sorted_nodup_rg.bam

### Cleanup
mv ${name}_nomulti_proper_sorted_nodup.txt $cpath
mv ${name}_nomulti_proper_sorted_nodup_rg.bam $cpath
mv ${name}_nomulti_proper_sorted_nodup_rg.bam.bai $cpath
rm -rf $temp_dir
