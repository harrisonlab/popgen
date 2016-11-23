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

for reads in $input/*/*/all_reads_trim_1.fastq.gz
do
    Jobs=$(qstat | grep 'sub_trinit' | wc -l)
    while [ $Jobs -gt 1 ]
    do
        output_dir=$(echo $reads | awk -F/ '{print $(NF-1)}')
        sleep 100
        printf "."
        Jobs=$(qstat | grep 'sub_trinit' | wc -l)
    done
reads2=$(echo "$reads" | sed 's/_1.fastq.gz/_2.fastq.gz/')
qsub $scripts/sub_trinity_assembly.sh $reads $reads2 $output_dir
done


#After done with that, carry out a tblastx BLAST comparison with the R genes
#in the old onion dataset to identify any potential new targets for bait design.

#If such are found, then use tblastn to check for the percentage identity in the sequences
#of R genes which already have homologs in both transcriptome datasets.
