#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 6
#$ -l virtual_free=1.5G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace12.blacklace

# Align raw reads to a pacbio assembly and then use this alignmeant to correct
# indels and substitutions in the assembly. Input: gzipped reads.

# ---------------
# Step 1
# Collect inputs
# ---------------

Assembly=$1
Read_F1=$2
Read_R1=$3
Read_F2=$4
Read_R2=$5
Read_F3=$6
Read_R3=$7

echo "Running Quickmerge with the following inputs:"
echo "Pacbio assembly - $Assembly"
echo "Forward trimmed reads - $Read_F1 $Read_F2 $Read_F3"
echo "Reverse trimmed reads - $Read_R1 $Read_R2 $Read_R3"

# ---------------
# Step 3
# Align seq reads
# ---------------
# Prepare the assembly for alignment
# Align reads against the assembly
# Convert the SAM file to BAM in preparation for sorting.
# Sort the BAM file, in preparation for SNP calling:
# Index the bam file

for i in $Read_F1 $Read_F2 $Read_F3 $Read_R1 $Read_R2 $Read_R3
do
file=$(basename $i)
gzip -dc < $i > $PWD/${file%.gz}
done

Read_R1_b=$(basename $Read_R1)
Read_R2_b=$(basename $Read_R2)
Read_R3_b=$(basename $Read_R3)

Read_F1_b=$(basename $Read_F1)
Read_F2_b=$(basename $Read_F2)
Read_F3_b=$(basename $Read_F3)

bowtie2-build $Assembly $Assembly.indexed
bowtie2 -p 6 -x $Assembly.indexed -1 ${Read_F1_b%.gz},${Read_F2_b%.gz},${Read_F3_b%.gz} -2 ${Read_R1_b%.gz},${Read_R2_b%.gz},${Read_R3_b%.gz} -S "$Assembly"_aligned.sam
samtools view -bS "$Assembly"_aligned.sam -o "$Assembly"_aligned.bam
samtools sort "$Assembly"_aligned.bam "$Assembly"_aligned_sorted
samtools index "$Assembly"_aligned_sorted.bam

# ---------------
# Step 4
# Run Pilon
# ---------------
# Run pilon to polish

JavaDir=/home/armita/prog/pilon
java -jar $JavaDir/pilon-1.17.jar --threads 6 --tracks --vcf --mindepth 30 --diploid --genome $Assembly  --frags "$Assembly"_aligned_sorted.bam --outdir $PWD/pilon

