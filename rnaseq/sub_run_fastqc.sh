#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l virtual_free=1G
#$ -l h=blacklace02.blacklace|blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

INFILE=$1
READS=$(basename $INFILE)

CUR_PATH=$PWD
WORK_DIR=$TMPDIR

mkdir -p $WORK_DIR
cp $INFILE $WORK_DIR    
cd $WORK_DIR

fastqc=/home/master_files/prog_master/bin/fastqc
$fastqc --nogroup $READS

cp -r $WORK_DIR/*zip $CUR_PATH

