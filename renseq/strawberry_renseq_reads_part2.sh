#!/bin/bash
scripts=/home/sobczm/bin/popgen/renseq
input=/home/sobczm/files/Helen_Bates_EMR.RH.ENQ-1704.A.01

#Carried out on blacklace11
###Trim Illumina adaptors from raw PacBio reads using a script from Giolai et al. (2016)
cd $input/Raw_reads_S1/F06_1/Analysis_Results
for a in *.bax.h5; do $scripts/trim-h5.py $a; done
cd $input/Raw_reads_S2/G06_1/Analysis_Results
for a in *.bax.h5; do $scripts/trim-h5.py $a; done

####Note: this analysis carried out on NIAB's triticum.
#Prepare CCS at 99% accuracy from minimum 3 reads pass
bb=/home/sobczm/bin/pitchfork/workspace/bax2bam/bin
ccs=/home/sobczm/bin/unanimity/build/ccs 
#Convert Bax to BAM for CCS
#S1
cd /home/sobczm/files/Helen_Bates_EMR.RH.ENQ-1704.A.01/Raw_reads_S1/F06_1/Analysis_Results
LD_LIBRARY_PATH=/home/sobczm/bin/pitchfork/deployment/lib $bb/bax2bam -o S1.bam m170127_051014_42165_c101119912550000001823252505221732_s1_p0.1.bax.h5 m170127_051014_42165_c101119912550000001823252505221732_s1_p0.2.bax.h5 m170127_051014_42165_c101119912550000001823252505221732_s1_p0.3.bax.h5
$ccs --minLength=1000 --minPredictedAccuracy=0.99 S1.bam.subreads.bam S1_ccs_3_99.bam
#S2
cd /home/sobczm/files/Helen_Bates_EMR.RH.ENQ-1704.A.01/Raw_reads_S2/G06_1/Analysis_Results
LD_LIBRARY_PATH=/home/sobczm/bin/pitchfork/deployment/lib $bb/bax2bam -o S2.bam m170127_093129_42165_c101119912550000001823252505221733_s1_p0.1.bax.h5 m170127_093129_42165_c101119912550000001823252505221733_s1_p0.2.bax.h5 m170127_093129_42165_c101119912550000001823252505221733_s1_p0.3.bax.h5
$ccs --minLength=1000 --minPredictedAccuracy=0.99 S2.bam.subreads.bam S2_ccs_3_99.bam
##CCS generated from only 20% sequences, in the majority of cases (50%) due to insufficient number of complete passes. Try to lower down to 2 passes, and see what happens.
$ccs --minLength=1000 --minPasses=2 --minPredictedAccuracy=0.99 S2.bam.subreads.bam S2_ccs_2_99.bam
#No difference observed so sticking with 3 passes.
#Copy the files over to the EMR cluster, convert to fastq and then fasta.
scp -r sobczm@10.1.10.170:/home/sobczm/files/Helen_Bates_EMR.RH.ENQ-1704.A.01/Raw_reads_S2/G06_1/Analysis_Results/S2_ccs_3_99.bam ./
scp -r sobczm@10.1.10.170:/home/sobczm/files/Helen_Bates_EMR.RH.ENQ-1704.A.01/Raw_reads_S1/F06_1/Analysis_Results/S1_ccs_3_99.bam ./

bedtools bamtofastq -i S1_ccs_3_99.bam -fq S1_ccs_3_99.fastq
bedtools bamtofastq -i S2_ccs_3_99.bam -fq S2_ccs_3_99.fastq
/home/sobczm/bin/seqtk/seqtk seq -a S1_ccs_3_99.fastq > S1_ccs_3_99.fasta
/home/sobczm/bin/seqtk/seqtk seq -a S2_ccs_3_99.fastq > S2_ccs_3_99.fasta
#Carry out minor variant calling with PacBio minorseq pipeline
#using canu and Falcon polished assembly as reference.

#Copy over HGAP assembly of reads
cd $input2/assembly
scp -r sobczm@10.1.10.170:/data/software/smrtanalysis/userdata/jobs/016/016448 ./
scp -r sobczm@10.1.10.170:/data/software/smrtanalysis/userdata/jobs/016/016449 ./

#Analyze contigs from HGAP
mv ./016448/data/polished_assembly.fasta.gz S1_HGAP_polished_assembly.fasta.gz
gunzip -d S1_HGAP_polished_assembly.fasta.gz

mv ./016449/data/polished_assembly.fasta.gz S2_HGAP_polished_assembly.fasta.gz
gunzip -d S2_HGAP_polished_assembly.fasta.gz

####STOP
####Failed attempts at other analyses follow below:

##Polish the RG Canu assembly with Pilon using Illumina reads,
# [Carried out on EMR cluster]
reads=/home/groups/harrisonlab/project_files/fragaria_x_ananassa/octoseq/PE
input2=/home/sobczm/popgen/renseq/strawberry/reads/Helen_Bates_EMR.RH.ENQ-1704.A.01
qsub $scripts/sub_pilon_3_libs.sh $input2/assembly/G06_assembly/G06.contigs.fasta $reads/1628_LIB19365_LDI16700_GCCAAT_L001_R1.fastq.gz $reads/1628_LIB19365_LDI16700_GCCAAT_L001_R2.fastq.gz $reads/1650_LIB19365_LDI16700_GCCAAT_L001_R1.fastq.gz $reads/1650_LIB19365_LDI16700_GCCAAT_L001_R2.fastq.gz $reads/1670_LIB19365_LDI16700_GCCAAT_L001_R1.fastq.gz $reads/1670_LIB19365_LDI16700_GCCAAT_L001_R2.fastq.gz 

#Get the following error from bowtie:
#Saw ASCII character 0 but expected 33-based Phred qual.
#terminate called after throwing an instance of 'int'
#Aborted (core dumped)
#(ERR): bowtie2-align exited with value 134
cd $input2/assembly/G06_assembly/
qsub $scripts/sub_pilon_3_libs.sh G06.contigs.fasta 1628_LIB19365_LDI16700_GCCAAT_L001_R1.fastq 1628_LIB19365_LDI16700_GCCAAT_L001_R2.fastq 1650_LIB19365_LDI16700_GCCAAT_L001_R1.fastq 1650_LIB19365_LDI16700_GCCAAT_L001_R2.fastq 1670_LIB19365_LDI16700_GCCAAT_L001_R1.fastq 1670_LIB19365_LDI16700_GCCAAT_L001_R2.fastq 
##But ASCII character 0 is @ which is within PHRED32.
##grepping for "M" or "d", one of the characters within PHRED64 brings no results. The file is clearly PHRED 32.

#Going to use bwa mem instead. Each library will be aligned separately and then merged with samtools merge.
qsub $scripts/sub_bwa_mem.sh Illumina 1628 G06.contigs.fasta 1628_LIB19365_LDI16700_GCCAAT_L001_R1.fastq 1628_LIB19365_LDI16700_GCCAAT_L001_R2.fastq 
#OK, now bwa mem complaining about not matching read names. Need to resort.
#sobczm@bio72:~/popgen/renseq/strawberry/reads/Helen_Bates_EMR.RH.ENQ-1704.A.01/assembly/G06_assembly$ grep -n D3NG5HQ1:391:H7HWCBCXX:1:1102:12271:51761 1628_LIB19365_LDI16700_GCCAAT_L001_R1.fastq 
#6779115:@D3NG5HQ1:391:H7HWCBCXX:1:1102:12271:51761 1:N:0:GCCAAT

#sobczm@bio72:~/popgen/renseq/strawberry/reads/Helen_Bates_EMR.RH.ENQ-1704.A.01/assembly/G06_assembly$ grep -n D3NG5HQ1:391:H7HWCBCXX:1:1102:12271:51761 1628_LIB19365_LDI16700_GCCAAT_L001_R2.fastq 
#10432589:@D3NG5HQ1:391:H7HWCBCXX:1:1102:12271:51761 2:N:0:GCCAAT

###TOO SLOW AND LOADS EVERYTHING INTO MEMORY:
qsub $scripts/sub_pairfq.sh 1628_LIB19365_LDI16700_GCCAAT_L001_R1.fastq 1628_LIB19365_LDI16700_GCCAAT_L001_R2.fastq

#Use FASTQ Tools instead
qsub $scripts/sub_fastq_sort.sh 1628_LIB19365_LDI16700_GCCAAT_L001_R1.fastq 

#Get an error involving too small bufffer size, even if mem allocated
#100 GB.
#Going to try something else:
qsub $scripts/sub_ngsutils_sort.sh 1628_LIB19365_LDI16700_GCCAAT_L001_R1.fastq
#Another error!!
#IOError: [Errno 24] Too many open files: '/home/groups/harrisonlab/project_files/popgen/renseq/strawberry/reads/Helen_Bates_EMR.RH.ENQ-1704.A.01/assembly/G06_assembly/.tmp.1628_LIB19365_LDI16700_GCCAAT_L001_R1.fastqVFHs1T'

#Sort using bash commands on blacklace11
cat 1628_LIB19365_LDI16700_GCCAAT_L001_R1.fastq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' >1628_LIB19365_LDI16700_GCCAAT_L001_R1_sort.fastq
cat 1628_LIB19365_LDI16700_GCCAAT_L001_R2.fastq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' >1628_LIB19365_LDI16700_GCCAAT_L001_R2_sort.fastq

#Still not sorted properly.

#Assemble each cultivar with Falcon on triticum.
#Load Falcon
falcon=/home/sobczm/bin/FALCON-integrate/
source /home/sobczm/bin/FALCON-integrate/env.sh

#Lambda test dataset
input=/data/projects/sobczm/lambda
cd $input
fc_run.py fc_run.cfg 

#Renseq F06
cd /data/projects/sobczm/renseq/F06
fc_run.py fc_run.cfg 
#Get the following error in /data/projects/sobczm/renseq/F06/0-rawreads/pwatcher.dir/stderr
#AssertionError: Not enough genome coverage (target=3000000000 < actual=1377800802)
###Changed genome size to 100000 and seed coverage to 10
