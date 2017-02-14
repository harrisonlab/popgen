#!/bin/bash
scripts=/home/sobczm/bin/popgen/renseq
input=/home/sobczm/popgen/renseq/strawberry/reads/Helen_Bates_EMR.RH.ENQ-1704.A.01

#Concatenate the raw reads
#Cultivar1
mkdir -p $input/assembly
cd $input/Raw_reads_S1/F06_1/Analysis_Results
for a in *subreads.fastqcd $input/assembly
docd $input/assembly
cat $a >> $input/assembly/F06_1_S1.fastq
done

#Cultivar2
cd $input/Raw_reads_S1/G06_1/Analysis_Results
for a in *subreads.fastq
do
cat $a >> $input/assembly/G06_1_S2.fastq
done

#Assemble each cultivar with Canu
cd $input/assembly
qsub $scripts/sub_canu_renseq.sh F06_1_S1.fastq 1m F06 F06_assembly
qsub $scripts/sub_canu_renseq.sh G06_1_S2.fastq 1m G06 G06_assembly

#Analyse HQ ROI, at 90% accuracy, min 1 full-length pass for NBS genes.
mkdir $input/assembly/nbs-parser
cp $input/Helen_Bates_EMR.RH.ENQ-1704.A.01_S1_MinFullPasses1_Accuracy90/data/reads_of_insert.fasta \
$input/assembly/nbs-parser/F06_1_S1_roi_fl1_a90.fasta
sh $scripts/sub_nlrparser.sh F06_1_S1_roi_fl1_a90.fasta
cp $input/Helen_Bates_EMR.RH.ENQ-1704.A.01_S2_MinFullPasses1_Accuracy90/data/reads_of_insert.fasta \
$input/assembly/nbs-parser/G06_1_S2_roi_fl1_a90.fasta
sh $scripts/sub_nlrparser.sh G06_1_S2_roi_fl1_a90.fasta

#Analyse the assemblies for NBS genes.
cd $input/assembly/F06_assembly
cat F06.contigs.fasta >> ../F06.assembly.fasta 
cat F06.unassembled.fasta >> ../F06.assembly.fasta 

cd $input/assembly/G06_assembly
cat G06.contigs.fasta >> ../G06.assembly.fasta 
cat G06.unassembled.fasta >> ../G06.assembly.fasta 
cp $input/assembly/*.assembly.fasta $input/assembly/nbs-parser

cd $input/assembly/nbs-parser
sh $scripts/sub_nlrparser.sh F06.assembly.fasta 
sh $scripts/sub_nlrparser.sh G06.assembly.fasta