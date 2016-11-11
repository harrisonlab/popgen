#$ -S /bin/bash
#$ -cwd
#$ -pe smp 20
#$ -l h_vmem=15.7G
#$ -l h=blacklace11.blacklace

left_reads=$1
right_reads=$2

#cpath=$PWD

#temp_dir="$TMPDIR"
#mkdir -p $temp_dir

#cp $left_reads $temp_dir
#cp $right_reads $temp_dir

lr=$(basename "$left_reads")
rr=$(basename "$right_reads")

#cd $temp_dir

trinity=/home/sobczm/bin/trinityrnaseq-2.2.0/Trinity
$trinity --seqType fa --max_memory 300G --CPU 20 --min_contig_length 300 --left $lr --right $rr

#rm *.fastq
#cp -r * $cpath
#rm -rf $temp_dir
