#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l virtual_free=1.0G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace11.blacklace

#Run read counts to be input into limma.
strandedness=$1
input_bam=$2
input_gff=$3
output_filepath=$4

htseq_count=/home/sobczm/.local/bin/htseq-count
$htseq_count --format=bam --stranded=$strandedness --order=pos --mode=union --type=exon --idattr=Parent $input_bam $input_gff > $output_filepath