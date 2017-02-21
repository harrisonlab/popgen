#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l virtual_free=15G
#$ -l h=blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

#Usage: qsub sub_pairfq.sh forward_reads reverse_reads

#The script accepts gunzipped (*.gz) FASTQ files: forward and reverse reads.
#It re-pairs paired-end reads (output: suffix _p) which are out of sync, and outputs two files with singletons (output: suffix _s), which are subsequently gunzipped (compressed) again.

forward=$1
reverse=$2
cpath=$PWD

temp_dir="$TMPDIR"
mkdir -p $temp_dir

cp -r $forward $reverse $temp_dir
input_f=$(basename "$forward")
input_r=$(basename "$reverse")
cd $temp_dir

gzip -d $input_f
gzip -d $input_r

unzipped_f=${$input_f%.gz}
unzipped_g=${$input_g%.gz}

pairfq=/home/sobczm/bin/Pairfq/bin/pairfq

perl $pairfq makepairs -f $unzipped_f -r $unzipped_g -fp ${unzipped_f%.f*}_p.fastq -rp ${unzipped_r%.f*}_p.fastq -fs ${unzipped_f%.f*}_s.fastq -rs ${unzipped_r%.f*}_s.fastq

gzip ${unzipped_f%.f*}_p.fastq 
gzip ${unzipped_r%.f*}_p.fastq 
gzip ${unzipped_f%.f*}_s.fastq
gzip ${unzipped_r%.f*}_s.fastq

cp -r *.gz $cpath
rm -rf $temp_dir