#!/bin/bash
input=/home/sobczm/popgen/renseq/input/reads
scripts=/home/sobczm/bin/popgen/renseq
fastqdump=/home/sobczm/bin/sratoolkit.2.8.0-ubuntu64/bin/fastq-dump

cd $input/common_onion

#Download all the reads in PHRED33 format.

#Read Download
#Common onion
#Han2016
$fastqdump --outdir han --gzip --skip-technical --readids --dumpbase  \
--split-files --clip SRR2814822
$fastqdump --outdir han --gzip --skip-technical --readids --dumpbase  \
--split-files --clip SRR2814821
$fastqdump --outdir han --gzip --skip-technical --readids --dumpbase  \
--split-files --clip SRR2814820
$fastqdump --outdir han --gzip --skip-technical --readids --dumpbase  \
--split-files --clip SRR2814819
$fastqdump --outdir han --gzip --skip-technical --readids --dumpbase  \
--split-files --clip SRR2814818
$fastqdump --outdir han --gzip --skip-technical --readids --dumpbase  \
--split-files --clip SRR2814817
$fastqdump --outdir han --gzip --skip-technical --readids --dumpbase  \
--split-files --clip SRR2814816
$fastqdump --outdir han --gzip --skip-technical --readids --dumpbase  \
--split-files --clip SRR2814815

#Cornell
$fastqdump --outdir cornell --gzip --skip-technical --readids --dumpbase  \
--split-files --clip SRR4418772
$fastqdump --outdir cornell --gzip --skip-technical --readids --dumpbase  \
--split-files --clip SRR4418771
$fastqdump --outdir cornell --gzip --skip-technical --readids --dumpbase  \
--split-files --clip SRR4418770
$fastqdump --outdir cornell --gzip --skip-technical --readids --dumpbase  \
--split-files --clip SRR4418769
$fastqdump --outdir cornell --gzip --skip-technical --readids --dumpbase  \
--split-files --clip SRR4418768
$fastqdump --outdir cornell --gzip --skip-technical --readids --dumpbase  \
--split-files --clip SRR4418767

#SP3B
$fastqdump --outdir SP3B --gzip --skip-technical --readids --dumpbase  \
--split-files --clip SRR1312067

#H6
$fastqdump --outdir H6 --gzip --skip-technical --readids --dumpbase  \
--split-files --clip SRR1312066

## Now, for the Welsh onion
cd $input/welsh_onion
#Liu2014
$fastqdump --outdir liu --gzip --skip-technical --readids --dumpbase  \
--split-files --clip SRR1609126
$fastqdump --outdir liu --gzip --skip-technical --readids --dumpbase  \
--split-files --clip SRR1609976
#Sun2016
$fastqdump --outdir sun --gzip --skip-technical --readids --dumpbase  \
--split-files --clip SRR1023632

#Read QC

for reads in $input/*/*/*_1.fastq.gz
do
    Jobs=$(qstat | grep 'sub_read_q' | wc -l)
    while [ $Jobs -gt 5 ]
    do
        sleep 10
        printf "."
        Jobs=$(qstat | grep 'sub_read_q' | wc -l)
    done
reads2=$(echo "$reads" | sed 's/_1.fastq.gz/_2.fastq.gz/')
qsub $scripts/sub_read_qc.sh $reads $reads2
done

#Remove the original downloaded files
for reads in $input/*/*/*_1.fastq.gz
do
    reads2=$(echo "$reads" | sed 's/_1.fastq.gz/_2.fastq.gz/')
    rm $reads $reads2
done
