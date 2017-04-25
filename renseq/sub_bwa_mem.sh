#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 12
#$ -l virtual_free=1G
#$ -l h=blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

bwa=/home/sobczm/bin/speedseq/src/bwa/bwa

read_type=$1
prefix=$2
assembly=$3
F_reads=$4
R_reads=$5

output=${F_reads%.fastq}.bam

# Align the data
$bwa index $assembly

if [ "$read_type" = "Illumina" ]
then
    $bwa mem -t 12 $assembly $F_reads $R_reads | samtools view -S -b - >$output
elif [ "$read_type" = "Pacbio" ]
then
    $bwa mem -t 12 -x pacbio $assembly $F_reads | samtools view -S -b - >$output
else
    echo "No read type given!"
fi

# Sort the alignment
samtools sort $output ${output%.bam}_sort

### Add group and sample name (prefix)
bamaddrg=/home/sobczm/bin/freebayes/bamaddrg/bamaddrg
output_rg=${output%.bam}_sort_rg.bam
$bamaddrg -b ${output%.bam}_sort.bam -s $prefix -r $prefix >$output_rg

#index
samtools index $output_rg


#Annotate the assembly with RGaugury