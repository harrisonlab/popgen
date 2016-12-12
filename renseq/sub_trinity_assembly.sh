#$ -S /bin/bash
#$ -cwd
#$ -pe smp 20
#$ -l h_vmem=15.7G
#$ -l h=blacklace11.blacklace

left_reads=$1
right_reads=$2
output_dir=$3

cpath=$PWD

temp_dir="$TMPDIR"
mkdir -p $temp_dir

cp $left_reads $temp_dir
cp $right_reads $temp_dir

lr=$(basename "$left_reads")
rr=$(basename "$right_reads")

cd $temp_dir

#First convert fastq to fasta.
fastool=/home/sobczm/bin/trinityrnaseq-2.2.0/trinity-plugins/fastool/fastool
zcat $lr | $fastool --illumina-trinity >all_reads_f.fasta
zcat $rr | $fastool --illumina-trinity >all_reads_r.fasta

trinity=/home/sobczm/bin/trinityrnaseq-2.2.0/Trinity
$trinity --seqType fa --max_memory 300G --CPU 20 --min_contig_length 300 \
--left all_reads_f.fasta --right all_reads_r.fasta


rm *.fastq
rm *.fastq.gz
rm all_reads_f.fasta
rm all_reads_r.fasta
cp -r * $cpath
mv $cpath/trinity_out_dir $cpath/$output_dir
rm -rf $temp_dir
