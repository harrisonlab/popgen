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

#Need to change the read ID headers to add /1 and /2 suffixes for downstream assembly
#in Trinity
for f in $input/*/*/*_1_trim*.gz
do
mv "$f" "$f~" &&
gzip -cd "$f~" | cut -d' ' -f1 | sed 's/\.1$/\/1/g' | gzip > "$f"
done

for f in $input/*/*/*_2_trim*.gz
do
mv "$f" "$f~" &&
gzip -cd "$f~" | cut -d' ' -f1 | sed 's/\.2$/\/2/g' | gzip > "$f"
done

#Remove temporary files
for f in $input/*/*/*~
do
rm $f
done


#Combine all trimmed reads for a given dataset into one forward and one reverse files
#Cornell
cd $input/common_onion/cornell
cat SRR4418767_1_trim.fq.gz SRR4418768_1_trim.fq.gz SRR4418769_1_trim.fq.gz \
SRR4418770_1_trim.fq.gz SRR4418771_1_trim.fq.gz SRR4418772_1_trim.fq.gz \
>all_reads_trim_1.fastq.gz
cat SRR4418767_2_trim.fq.gz SRR4418768_2_trim.fq.gz SRR4418769_2_trim.fq.gz \
SRR4418770_2_trim.fq.gz SRR4418771_2_trim.fq.gz SRR4418772_2_trim.fq.gz \
>all_reads_trim_2.fastq.gz
#Han2016
cd $input/common_onion/han
cat SRR2814815_1_trim.fq.gz SRR2814822_1_trim.fq.gz >all_reads_trim_1.fastq.gz
cat SRR2814815_2_trim.fq.gz SRR2814822_2_trim.fq.gz >all_reads_trim_2.fastq.gz
cd $input/common_onion/H6
mv SRR1312066_1_trim.fq.gz all_reads_trim_1.fastq.gz
mv SRR1312066_2_trim.fq.gz all_reads_trim_2.fastq.gz
cd $input/common_onion/SP3B
mv SRR1312067_1_trim.fq.gz all_reads_trim_1.fastq.gz
mv SRR1312067_2_trim.fq.gz all_reads_trim_2.fastq.gz
cd $input/welsh_onion/liu
cat SRR1609126_1_trim.fq.gz SRR1609976_1_trim.fq.gz >all_reads_trim_1.fastq.gz
cat SRR1609126_2_trim.fq.gz SRR1609976_2_trim.fq.gz >all_reads_trim_2.fastq.gz
cd $input/welsh_onion/sun
mv SRR1023632_1_trim.fq.gz all_reads_trim_1.fastq.gz
mv SRR1023632_2_trim.fq.gz all_reads_trim_2.fastq.gz

#Remove the original downloaded files
for reads in $input/*/*/*_1.fastq.gz
do
    reads2=$(echo "$reads" | sed 's/_1.fastq.gz/_2.fastq.gz/')
    rm $reads $reads2
done
