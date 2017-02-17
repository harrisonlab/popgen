#!/bin/bash
input_hap=/home/groups/harrisonlab/project_files/venturia/qc_dna/paired/v.inaequalis
input_hap_assembly=/home/groups/harrisonlab/project_files/venturia/repeat_masked/v.inaequalis/172_pacbio/filtered_contigs_repmask/172_pacbio_contigs_unmasked.fa

input_dip=/home/groups/harrisonlab/project_files/phytophthora_fragariae/qc_dna/paired/P.fragariae
input_dip_assembly=/home/groups/harrisonlab/project_files/phytophthora_fragariae/repeat_masked/P.fragariae/Bc16/filtered_contigs_repmask/95m_contigs_unmask

scripts=/home/sobczm/bin/popgen/snp
lumpy=/home/sobczm/bin/lumpy-sv/bin

toyset=/home/sobczm/bin/speedseq/example/toyset
cd $toyset

qsub $scripts/sub_bwa_mem.sh human_g1k_v37_20_42220611-42542245.fasta NA12878.20slice.30X_1.fastq NA12878.20slice.30X_2.fastq Pacbio test

qsub $scripts/sub_lumpy.sh NA12878.20slice.30X_1_rg.bam NA12878.20slice.30X_1_rg_splitters.bam NA12878.20slice.30X_1_rg_discordants.bam

