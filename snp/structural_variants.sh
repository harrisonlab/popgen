#!/bin/bash
#Structural variants (duplications, deletions, inversions, translocations) detection using read-pair configuration, split-reads, read-depth.
input_hap=/home/groups/harrisonlab/project_files/venturia/qc_dna/paired/v.inaequalis
input_hap_assembly=/home/groups/harrisonlab/project_files/venturia/repeat_masked/v.inaequalis/172_pacbio/filtered_contigs_repmask/172_pacbio_contigs_unmasked.fa

input_dip=/home/groups/harrisonlab/project_files/phytophthora_fragariae/qc_dna/paired/P.fragariae
input_dip_assembly=/home/groups/harrisonlab/project_files/phytophthora_fragariae/repeat_masked/P.fragariae/Bc16/filtered_contigs_repmask/95m_contigs_unmasked.fa

scripts=/home/sobczm/bin/popgen/snp
input=/home/sobczm/popgen/snp/sv_calling

#Toyset example
toyset=/home/sobczm/bin/speedseq/example/toyset
cd $toyset

#Read alignment with bwa-mem
qsub $scripts/sub_bwa_mem.sh human_g1k_v37_20_42220611-42542245.fasta NA12878.20slice.30X_1.fastq NA12878.20slice.30X_2.fastq Pacbio test
#Structural variant calling
qsub $scripts/sub_lumpy.sh NA12878.20slice.30X_1_rg.bam NA12878.20slice.30X_1_rg_splitters.bam NA12878.20slice.30X_1_rg_discordants.bam

mkdir -p $input/pfrag
cd $input/pfrag
###Phytophthora fragariae (diploid) example
#Concatenated PaCBio subreads from one sample
pacbio_con_1=/home/groups/harrisonlab/project_files/phytophthora_fragariae/raw_dna/pacbio/P.fragariae/Bc16/extracted
qsub $scripts/sub_bwa_mem.sh Pacbio vinequalis $input_dip_assembly $pacbio_con_1/concatenated_pacbio.fastq 

#Illumina samples
#QC-trimmed reads. Compressed here but don't have to be
for sample in $input_dip/*
do
reads_forward=$sample/F/*trim.fq.gz
reads_reverse=$sample/R/*trim.fq.gz
qsub $scripts/sub_bwa_mem.sh Illumina pfrag $input_dip_assembly $reads_forward $reads_reverse 
done


mkdir -p $input/vinequalis
cd $input/vinequalis
###Venturia inequalis (haploid) example
#Illumina samples.
#QC-trimmed reads. Compressed here but don't have to be
for sample in $input_hap/*
do
reads_forward=$sample/F/*trim.fq.gz
reads_reverse=$sample/R/*trim.fq.gz
qsub $scripts/sub_bwa_mem.sh Illumina vinequalis $input_hap_assembly $reads_forward $reads_reverse 
done

#Concatenated PaCBio subreads from one sample
pacbio_con_2=/home/groups/harrisonlab/project_files/venturia/raw_dna/pacbio/v.inaequalis/172_pacbio/extracted
qsub $scripts/sub_bwa_mem.sh Pacbio vinequalis $input_hap_assembly $pacbio_con_2/concatenated_pacbio.fastq 
done