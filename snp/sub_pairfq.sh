#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l virtual_free=15G
#$ -l h=blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

#Usage: qsub sub_pairfq.sh forward_reads reverse_reads

#The script accepts FASTQ files: forward and reverse reads.
#It re-pairs paired-end reads (output suffix: _f_p.fastq and _r_p.fastq) which are out of sync, and also outputs two files with singletons (output suffix: _f_s.fastq and _r_s.fastq), which are subsequently gunzipped (compressed).

forward=$1
reverse=$2
cpath=$PWD

temp_dir="$TMPDIR"
mkdir -p $temp_dir

cp -r $forward $reverse $temp_dir
input_f=$(basename "$forward")
input_r=$(basename "$reverse")
cd $temp_dir

pairfq=/home/sobczm/bin/Pairfq/bin/pairfq

mv $input_f input_f.fastq
mv $input_r input_r.fastq

perl $pairfq makepairs -c gzip -f input_f.fastq -r input_r.fastq -fp f_p.fastq -rp r_p.fastq -fs f_s.fastq -rs r_s.fastq

mv f_p.fastq.gz "${input_f%.f*}_f_p.fastq.gz"
mv r_p.fastq.gz "${input_r%.f*}_r_p.fastq.gz"
mv f_s.fastq.gz "${input_f%.f*}_f_s.fastq.gz"
mv r_s.fastq.gz "${input_r%.f*}_r_s.fastq.gz"

cp -r *.gz $cpath
rm -rf $temp_dir