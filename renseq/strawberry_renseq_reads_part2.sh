#!/bin/bash
scripts=/home/sobczm/bin/popgen/renseq
input=/home/sobczm/files/Helen_Bates_EMR.RH.ENQ-1704.A.01

####Note: analysis carried out on NIAB's triticum.

###Trim Illumina adaptors from raw PacBio reads using a script from Giolai et al. (2016)
cd $input/Raw_reads_S1/F06_1/Analysis_Results
for a in *.bax.h5; do $scripts/trim-h5.py $a; done
cd $input/Raw_reads_S2/G06_1/Analysis_Results
for a in *.bax.h5; do $scripts/trim-h5.py $a; done

##Polish the RG Canu assembly with Pilon using Illumina reads,
# [Carried out on EMR cluster]
reads=/home/groups/harrisonlab/project_files/fragaria_x_ananassa/octoseq/PE
input2=/home/sobczm/popgen/renseq/strawberry/reads/Helen_Bates_EMR.RH.ENQ-1704.A.01
qsub $scripts/sub_pilon_3_libs.sh $input2/assembly/G06_assembly/G06.contigs.fasta $reads/1628_LIB19365_LDI16700_GCCAAT_L001_R1.fastq.gz $reads/1628_LIB19365_LDI16700_GCCAAT_L001_R2.fastq.gz $reads/1650_LIB19365_LDI16700_GCCAAT_L001_R1.fastq.gz $reads/1650_LIB19365_LDI16700_GCCAAT_L001_R2.fastq.gz $reads/1670_LIB19365_LDI16700_GCCAAT_L001_R1.fastq.gz $reads/1670_LIB19365_LDI16700_GCCAAT_L001_R2.fastq.gz 

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
#For RedGauntlent, polish with Illumina reads using Canu.


#Carry out minor variant calling with PacBio minorseq pipeline
#using canu and Falcon polished assembly as reference.
