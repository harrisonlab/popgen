#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 16
#$ -l virtual_free=1G
#$ -l h=blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

bwa=/home/sobczm/bin/speedseq/src/bwa/bwa
scripts=/home/sobczm/bin/lumpy-sv/scripts

assembly=$1
F_reads=$2
R_reads=$3
read_type=$4
prefix=$5

cpath=$PWD

temp_dir="$TMPDIR"
mkdir -p $temp_dir

cp -r $F_reads $R_reads $assembly $temp_dir
fr=$(basename "$F_reads")
rr=$(basename "$R_reads")
as=$(basename "$assembly")

cd $temp_dir
output=${fr%.*}.bam

# Align the data
$bwa index $as

if [ "$read_type" = "Illumina" ]
then
    $bwa mem -t 16 $as $fr $rr | samtools view -S -b - >$output
elif [ "$read_type" = "Pacbio" ]
then
    $bwa mem -t 16 -x pacbio $as $fr $rr | samtools view -S -b - >$output
else
    echo "No read type given!"
fi


picard=/home/sobczm/bin/picard-tools-2.5.0/picard.jar
output_rg=${output%.bam}_rg.bam
### Add group and sample name (prefix)
java -jar $picard AddOrReplaceReadGroups INPUT=$output OUTPUT=$output_rg SORT_ORDER=coordinate CREATE_INDEX=true RGID=$prefix  RGSM=$prefix  RGPL=Illumina RGLB=library RGPU=barcode VALIDATION_STRINGENCY=LENIENT 
samtools index $output_rg 


# Extract the discordant paired-end alignments.
samtools view -b -F 1294 $output_rg > ${output_rg%.bam}_discordants_unsorted.bam

# Extract the split-read alignments
samtools view -h $output_rg | $scripts/extractSplitReads_BwaMem -i stdin | samtools view -Sb - > ${output_rg%.bam}_splitters_unsorted.bam

# Sort both alignments
samtools sort ${output_rg%.bam}_discordants_unsorted.bam ${output_rg%.bam}_discordants
samtools sort ${output_rg%.bam}_splitters_unsorted.bam ${output_rg%.bam}_splitters

rm ${output_rg%.bam}_splitters_unsorted.bam ${output_rg%.bam}_discordants_unsorted.bam $output
rm $as $fr $rr
cp -r * $cpath
rm -rf $temp_dir