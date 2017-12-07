###Copy the input files to triticum
datadir=/data/seq_data/external/20171110_EmilyFenellaRenseq_pacbio

##Carried out on blacklace11
###Trim Illumina adaptors from raw PacBio reads using a script from Giolai et al. (2016)
input=/home/sobczm/popgen/renseq/strawberry/reads/pacbio_emxfe
scripts=/home/sobczm/bin/popgen/renseq

tar -zvxf $datadir/Helen_Bates_EMR_RH_ENQ-1704_A_01_S3.tar.gz
tar -zxvf $datadir/Helen_Bates_EMR_RH_ENQ-1704_A_01_S4.tar.gz 

cd $input/Helen_Bates_EMR_RH_ENQ-1704_A_01_S3/2017_10_25_PSEQ1569_406/SAM32001_PRO1514_S3_HMWDNA_Emily/raw_reads/D01_1/Analysis_Results
for a in *.bax.h5; do $scripts/trim-h5.py $a; done

cd $input/Helen_Bates_EMR_RH_ENQ-1704_A_01_S4/2017_10_25_PSEQ1569_406/SAM32002_PRO1514_S4_HMWDNA_Fenella/raw_reads/E01_1/Analysis_Results
for a in *.bax.h5; do $scripts/trim-h5.py $a; done

scp -r $input/Helen_Bates_EMR_RH_ENQ-1704_A_01_S3 sobczm@10.1.10.170:/data/projects/sobczm/pacbio_emxfe
scp -r $input/Helen_Bates_EMR_RH_ENQ-1704_A_01_S4 sobczm@10.1.10.170:/data/projects/sobczm/pacbio_emxfe

#Prepare CCS at 99% accuracy from minimum 3 reads pass
bb=/home/sobczm/bin/pitchfork/workspace/bax2bam/bin
ccs=/home/sobczm/bin/unanimity/build/ccs 

#S3
cd /data/projects/sobczm/pacbio_emxfe/Helen_Bates_EMR_RH_ENQ-1704_A_01_S3/Analysis_Results
LD_LIBRARY_PATH=/home/sobczm/bin/pitchfork/deployment/lib $bb/bax2bam -o S3.bam m171026_073505_42165_c101248002550000001823283011021703_s1_p0.1.bax.h5 m171026_073505_42165_c101248002550000001823283011021703_s1_p0.2.bax.h5 m171026_073505_42165_c101248002550000001823283011021703_s1_p0.3.bax.h5
$ccs --minLength=1000 --minPredictedAccuracy=0.99 S3.bam.subreads.bam S3_ccs_3_99.bam

cd /data/projects/sobczm/pacbio_emxfe/Helen_Bates_EMR_RH_ENQ-1704_A_01_S4/Analysis_Results
LD_LIBRARY_PATH=/home/sobczm/bin/pitchfork/deployment/lib $bb/bax2bam -o S4.bam m171026_120100_42165_c101248002550000001823283011021704_s1_p0.1.bax.h5 m171026_120100_42165_c101248002550000001823283011021704_s1_p0.2.bax.h5 m171026_120100_42165_c101248002550000001823283011021704_s1_p0.3.bax.h5
$ccs --minLength=1000 --minPredictedAccuracy=0.99 S4.bam.subreads.bam S4_ccs_3_99.bam

#Copy the files over to the EMR cluster, convert to fastq and then fasta.
cd $input/analysis
scp -r sobczm@10.1.10.170:/data/projects/sobczm/pacbio_emxfe/Helen_Bates_EMR_RH_ENQ-1704_A_01_S3/Analysis_Results/S3_ccs_3_99.bam ./
scp -r sobczm@10.1.10.170:/data/projects/sobczm/pacbio_emxfe/Helen_Bates_EMR_RH_ENQ-1704_A_01_S4/Analysis_Results/S4_ccs_3_99.bam ./

bedtools bamtofastq -i S3_ccs_3_99.bam -fq S3_ccs_3_99.fastq
bedtools bamtofastq -i S4_ccs_3_99.bam -fq S4_ccs_3_99.fastq
/home/sobczm/bin/seqtk/seqtk seq -a S3_ccs_3_99.fastq > S3_ccs_3_99.fasta
/home/sobczm/bin/seqtk/seqtk seq -a S4_ccs_3_99.fastq > S4_ccs_3_99.fasta