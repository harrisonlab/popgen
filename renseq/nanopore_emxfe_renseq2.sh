#!/bin/bash
scripts=/home/sobczm/bin/popgen/renseq
input=/home/sobczm/popgen/renseq/strawberry/reads/albacore_emxfe
cd $input

##Second part of the analysis of Renseq data - following read assembly and correction. 
#Steps undertaken this far:
#Pipeline A
#1) Split reads and trim adapters using porechop
#2) Error correction with LoRMA -> result in polished HQ reads.

#Pipeline B 
#1) Split reads and trim adapters using porechop
#2) Read error correction with Canu
#3) Assembly of reads with smartdenovo
#4) 10 rounds of assembly polishing  using racon

#Current goals.
#Analyze LoRMA-corrected reads on their own.
#Call variants on the assembly from Pipeline B - using LoRMA-corrected reads as well as Canu-corrected reads as input.
#Illumina reads on hold now - not available yet. 

#### 1) Variant calling with GATK.
# - Align to reference with bwa-mem.
# - Add read groups
# - Don't mark duplicates as few reads available and PCR duplicates not likely. 
# - variant calling with GATK


###Map Canu-corrected reads 
qsub $scripts/sub_bwa_nanopolish.sh smartdenovo_barcode12_fenella_trimmed_all/racon/smartdenovo_barcode12_fenella_trimmed_all_racon_round_10.fasta barcode12_fenella_trimmed_all/barcode12_fenella_trimmed_all.trimmedReads.fasta.gz bwa_canu_fenella_all_own
qsub $scripts/sub_bwa_nanopolish.sh smartdenovo_barcode11_emily_trimmed_all/racon/smartdenovo_barcode11_emily_trimmed_all_racon_round_10.fasta barcode12_fenella_trimmed_all/barcode12_fenella_trimmed_all.trimmedReads.fasta.gz bwa_canu_fenella_all_cross
qsub $scripts/sub_bwa_nanopolish.sh smartdenovo_barcode11_emily_trimmed_all/racon/smartdenovo_barcode11_emily_trimmed_all_racon_round_10.fasta barcode11_emily_trimmed_all/barcode11_emily_trimmed_all.trimmedReads.fasta.gz bwa_canu_emily_all_own
qsub $scripts/sub_bwa_nanopolish.sh smartdenovo_barcode12_fenella_trimmed_all/racon/smartdenovo_barcode12_fenella_trimmed_all_racon_round_10.fasta barcode11_emily_trimmed_all/barcode11_emily_trimmed_all.trimmedReads.fasta.gz bwa_canu_emily_all_cross

###Map LoRMA-corrected reads 
qsub $scripts/sub_bwa_nanopolish.sh smartdenovo_barcode11_emily_trimmed_all/racon/smartdenovo_barcode11_emily_trimmed_all_racon_round_10.fasta lorma_barcode11_emily_all.fasta.gz bwa_lorma_emily_all_own
qsub $scripts/sub_bwa_nanopolish.sh smartdenovo_barcode12_fenella_trimmed_all/racon/smartdenovo_barcode12_fenella_trimmed_all_racon_round_10.fasta lorma_barcode11_emily_all.fasta.gz bwa_lorma_emily_all_cross
qsub $scripts/sub_bwa_nanopolish.sh smartdenovo_barcode12_fenella_trimmed_all/racon/smartdenovo_barcode12_fenella_trimmed_all_racon_round_10.fasta lorma_barcode12_fenella_all.fasta.gz bwa_lorma_fenella_all_own
qsub $scripts/sub_bwa_nanopolish.sh smartdenovo_barcode11_emily_trimmed_all/racon/smartdenovo_barcode11_emily_trimmed_all_racon_round_10.fasta lorma_barcode12_fenella_all.fasta.gz bwa_lorma_fenella_all_cross

#Prepare mappings for variant calling with GATK.
qsub $scripts/renseq_pre_snp_calling_prep.sh bwa_canu_fenella_all_own/reads.sorted.bam fenella $input/bwa_canu_fenella_all_own
qsub $scripts/renseq_pre_snp_calling_prep.sh bwa_canu_fenella_all_cross/reads.sorted.bam fenella $input/bwa_canu_fenella_all_cross
qsub $scripts/renseq_pre_snp_calling_prep.sh bwa_lorma_fenella_all_own/reads.sorted.bam fenella $input/bwa_lorma_fenella_all_own
qsub $scripts/renseq_pre_snp_calling_prep.sh bwa_lorma_fenella_all_cross/reads.sorted.bam fenella $input/bwa_lorma_fenella_all_cross
qsub $scripts/renseq_pre_snp_calling_prep.sh bwa_canu_emily_all_own/reads.sorted.bam emily $input/bwa_canu_emily_all_own
qsub $scripts/renseq_pre_snp_calling_prep.sh bwa_canu_emily_all_cross/reads.sorted.bam emily $input/bwa_canu_emily_all_cross
qsub $scripts/renseq_pre_snp_calling_prep.sh bwa_lorma_emily_all_own/reads.sorted.bam emily $input/bwa_lorma_emily_all_own
qsub $scripts/renseq_pre_snp_calling_prep.sh bwa_lorma_emily_all_cross/reads.sorted.bam emily $input/bwa_lorma_emily_all_cross

#Copy the assemblies over and prepare genome reference indexes required by GATK
cp smartdenovo_barcode11_emily_trimmed_all/racon/smartdenovo_barcode11_emily_trimmed_all_racon_round_10.fasta $input/gatk
cp smartdenovo_barcode12_fenella_trimmed_all/racon/smartdenovo_barcode12_fenella_trimmed_all_racon_round_10.fasta $input/gatk
java -jar /home/sobczm/bin/picard-tools-2.5.0/picard.jar CreateSequenceDictionary R=gatk/smartdenovo_barcode12_fenella_trimmed_all_racon_round_10.fasta O=gatk/smartdenovo_barcode12_fenella_trimmed_all_racon_round_10.dict
samtools faidx gatk/smartdenovo_barcode12_fenella_trimmed_all_racon_round_10.fasta
java -jar /home/sobczm/bin/picard-tools-2.5.0/picard.jar CreateSequenceDictionary R=gatk/smartdenovo_barcode11_emily_trimmed_all_racon_round_10.fasta O=gatk/smartdenovo_barcode11_emily_trimmed_all_racon_round_10.dict
samtools faidx gatk/smartdenovo_barcode11_emily_trimmed_all_racon_round_10.fasta

#Run GATK variant calling.
qsub $scripts/sub_gatk_canu_emily.sh
qsub $scripts/sub_gatk_canu_fenella.sh
qsub $scripts/sub_gatk_lorma_emily.sh
qsub $scripts/sub_gatk_lorma_fenella.sh