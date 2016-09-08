#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace11.blacklace|blacklace12.blacklace

assembly=$1
db=$2
blast=/home/sobczm/bin/popgen/codon/run_blast_reciprocal.pl

cpath=$PWD
temp_dir="$TMPDIR"

mkdir -p $temp_dir
cp $assembly $db* $temp_dir
cd $temp_dir

$blast $assembly $db

cp $temp_dir/*vs* $cpath
rm -r $temp_dir
