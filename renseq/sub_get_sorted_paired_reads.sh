#$ -S /bin/bash
#$ -cwd
#$ -pe smp 4
#$ -l h_vmem=4G
#$ -l h=blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

unmapped_bam=$1
mapped_bam=$2

cpath=$PWD

ma=$(basename "$mapped_bam")
un=$(basename "$unmapped_bam")

temp_dir="$TMPDIR"
mkdir -p $temp_dir
cp $mapped_bam $temp_dir
cp $unmapped_bam $temp_dir
cd $temp_dir

bamtools convert -in $un -out unmapped.fastq -format fastq
bamtools convert -in $ma -out mapped.fastq -format fastq

scripts=/home/sobczm/bin/popgen/renseq

python $scripts/find_matching_reads.py unmapped.fastq mapped.fastq f0.fastq r0.fastq
python $scripts/find_matching_reads2.py f0.fastq r0.fastq r.fastq
python $scripts/find_matching_reads2.py r.fastq f0.fastq f.fastq

cat r.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > r_sorted.fastq
cat f.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > f_sorted.fastq

cp -r *.fastq $cpath
#cp -r r_sorted.fastq $cpath
#cp -r f_sorted.fastq $cpath
rm -rf $temp_dir
