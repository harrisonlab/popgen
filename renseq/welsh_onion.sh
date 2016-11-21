#!/bin/bash
input=/home/sobczm/popgen/renseq/input/reads
scripts=/home/sobczm/bin/popgen/renseq
fastqdump=/home/sobczm/bin/sratoolkit.2.8.0-ubuntu64/bin/fastq-dump

#As the number of R genes discovered in the onion transcriptomes was quite low
#mining the other Allium species transcriptomes for them. Towards this end
#will download the reads for published Welsh onion (and also remaining common onion)
#transcriptomes and re-assamble them with Trinity.
cd $input
sh ./download_reads.sh
#Combine all trimmed reads for a given dataset into one forward and one reverse files
#Cornell
cd $input/common_onion/cornell
cat SRR4418767_1_trim.fq.gz SRR4418768_1_trim.fastq.gz SRR4418769_1_trim.fastq.gz \
SRR4418770_1_trim.fastq.gz SRR4418771_1_trim.fastq.gz SRR4418772_1_trim.fastq.gz \
>all_reads_trim_1.fastq.gz
cat SRR4418767_2_trim.fq.gz SRR4418768_2_trim.fastq.gz SRR4418769_2_trim.fastq.gz \
SRR4418770_2_trim.fastq.gz SRR4418771_2_trim.fastq.gz SRR4418772_2_trim.fastq.gz \
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

for reads in $input/*/*/all_reads_trim_1.fastq.gz
do
    Jobs=$(qstat | grep 'sub_trinit' | wc -l)
    while [ $Jobs -gt 1 ]
    do
        sleep 100
        printf "."
        Jobs=$(qstat | grep 'sub_trinit' | wc -l)
    done
reads2=$(echo "$reads" | sed 's/_1.fastq.gz/_2.fastq.gz/')
qsub $scripts/sub_trinity_assembly.sh $reads $reads2
done

sub_trinit
qsub


#After done with that, carry out a tblastx BLAST comparison with the R genes
#in the old onion dataset to identify any potential new targets for bait design.

#If such are found, then use tblastn to check for the percentage identity in the sequences
#of R genes which already have homologs in both transcriptome datasets.
