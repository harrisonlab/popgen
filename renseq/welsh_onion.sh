#!/bin/bash
input=/home/sobczm/popgen/renseq/input/reads/
scripts=/home/sobczm/bin/popgen/renseq
fastqdump=/home/sobczm/bin/sratoolkit.2.8.0-ubuntu64/bin/fastq-dump

#As the number of R genes discovered in the onion transcriptomes was quite low
#mining the other Allium species transcriptomes for them. Towards this end
#will download the reads for published Welsh onion (and also remaining common onion)
#transcriptomes and re-assamble them with Trinity.

cd $input/common_onion
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




#Read Assembly




#After done with that, carry out a tblastx BLAST comparison with the R genes
#in the old onion dataset to identify any potential new targets for bait design.

#If such are found, then use tblastn to check for the percentage identity in the sequences
#of R genes which already have homologs in both transcriptome datasets.
