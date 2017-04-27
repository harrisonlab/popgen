#$ -S /bin/bash
#$ -cwd
#$ -pe smp 4
#$ -l h_vmem=1G
#$ -l h=blacklace06.blacklace
###

#Input: Just the name for output results files. The script needs to be run from the directory containing all the output files from sub_bwa_mem.sh for the samples of interest (but not any others!!)

#Output: Genotype (*.gt) and VCF (*.vcf) files with structural variants.

output=$1

#First, obtain command line arguments with: 1) all the alignments
# 2) discordant alignments, 3) split alignments to be used as input
# by Lumpy

#1)
bams=()
for b in *rg.bam; do bams+=("$b"); done;
bam=$(IFS=, ; echo "${bams[*]}")

#2)
discs=()
for b in *discordants.bam; do discs+=("$b"); done;
discordant=$(IFS=, ; echo "${discs[*]}")

#3)
splits=()
for b in *splitters.bam; do splits+=("$b"); done;
splitters=$(IFS=, ; echo "${splits[*]}")

lumpy=/home/sobczm/bin/lumpy-sv/bin
out=${output}.gt

#Detect structural variants
$lumpy/lumpyexpress -B $bam -S $splitters -D $discordant -o $out

#Call genotypes
svtyper=/home/sobczm/bin/svtyper
$svtyper/svtyper -B $bam -S $splitters -i $out >${output}.vcf

