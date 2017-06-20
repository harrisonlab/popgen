#!/bin/bash
input=/home/sobczm/popgen/rnaseq/ananassa
scripts=/home/sobczm/bin/popgen/rnaseq

#Quality trim the raw sequence data.
Raw_Reads=/home/groups/harrisonlab/project_files/fragaria_x_ananassa/octoseq/RNAseq

for reads in $(ls $Raw_Reads/*1.fastq.gz); 
do
    Jobs=$(qstat | grep 'sub_read_q' | wc -l)
    while [ $Jobs -gt 10 ]
    do
        sleep 10
        printf "."
        Jobs=$(qstat | grep 'sub_read_q' | wc -l)
    done
reads2=$(echo "$reads" | sed 's/1.fastq.gz/2.fastq.gz/')
qsub $scripts/sub_read_qc.sh $reads $reads2
done