#$ -S /bin/bash
#$ -cwd
#$ -pe smp 8
#$ -l virtual_free=1.2G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace11.blacklace


Usage="sub_star.sh InGenome.fa InReadF.fa InReadR.fa Output_directory [GeneLocations.gff]"
echo $Usage

# ---------------
# Step 1
# Collect inputs
# ---------------

InGenome=$(basename $1)
echo $InGenome
InReadF=$(basename $2)
echo $InReadF
InReadR=$(basename $3)
echo $InReadR
out_dir=$4
echo $out_dir

# Set working directory
CurDir=$PWD
WorkDir=$TMPDIR/star
mkdir -p $WorkDir
mkdir -p $out_dir

# Copy over input files
cd $WorkDir
echo $PWD

cp -r $1 ./
cp $2 ./
cp $3 ./


# ---------------
# Step 2=3
# Run STAR
# ---------------

echo "Aligning RNAseq reads"

STAR \
--genomeDir $InGenome \
--readFilesCommand zcat \
--readFilesIn $InReadF $InReadR \
--outSAMtype BAM SortedByCoordinate \
--outFilterMultimapNmax 1 \
--alignIntronMin 20 \
--twopassMode Basic \
--outReadsUnmapped Fastx \
--runThreadN 8

samtools index Aligned.sortedByCoord.out.bam

rm -rf $WorkDir/$InGenome
rm $WorkDir/$InReadF
rm $WorkDir/$InReadR
cp -r $WorkDir/* $out_dir