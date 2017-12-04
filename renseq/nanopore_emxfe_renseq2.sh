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
#Analyze LoRMA-corrected reads on their own. (subsequently dropped due to low quality)
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

###Map LoRMA-corrected reads (subsequently dropped due to low quality)
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

#Run GATK variant calling (ploidy = 8).
qsub $scripts/sub_gatk_canu_emily.sh
qsub $scripts/sub_gatk_canu_fenella.sh
qsub $scripts/sub_gatk_lorma_emily.sh
qsub $scripts/sub_gatk_lorma_fenella.sh

#Run GATK variant calling (ploidy = 2).
qsub $scripts/sub_gatk_canu_emily2.sh
qsub $scripts/sub_gatk_canu_fenella2.sh
qsub $scripts/sub_gatk_lorma_emily2.sh
qsub $scripts/sub_gatk_lorma_fenella2.sh

#Input files to be analyzed.
cd $input/analysis
#Canu-corrected reads
#LoRMA corrected reads
#SMRT de novo assembly

#Blast the raw reads to establish hits R genes.
#Copy the blast db
cp -r /home/sobczm/popgen/renseq/strawberry/reads/Helen_Bates_EMR.RH.ENQ-1704.A.01/analysis/vesca_v1.1_nblrrs_augustus_mrna_nucl.db* ./

for a in barcode11_emily_trimmed_all.trimmedReads.fasta barcode12_fenella_trimmed_all.trimmedReads.fasta lorma_barcode11_emily_all.fasta lorma_barcode12_fenella_all.fasta smartdenovo_barcode11_emily_trimmed_all_racon_round_10.fasta smartdenovo_barcode12_fenella_trimmed_all_racon_round_10.fasta
do 
#Remove spaces in FASTA header, and substitute commas with _, as required by BLAST suite.
sed 's, ,_,g' -i $a
sed 's/,/_/g' -i $a
perl /home/sobczm/bin/popgen/other/assemblathon_stats.pl $a > ${a%.fasta}.stat
qsub $scripts/sub_blastn_renseq.sh $a vesca_v1.1_nblrrs_augustus_cds_nucl.db
sh $scripts/sub_nlrparser.sh $(basename $a)
done

#Create blast databases out of the assemblies
for assembly in barcode11_emily_trimmed_all.trimmedReads.fasta barcode12_fenella_trimmed_all.trimmedReads.fasta lorma_barcode11_emily_all.fasta lorma_barcode12_fenella_all.fasta smartdenovo_barcode11_emily_trimmed_all_racon_round_10.fasta smartdenovo_barcode12_fenella_trimmed_all_racon_round_10.fasta
do
makeblastdb -in $assembly -input_type fasta -dbtype nucl \
-title "${assembly%.*}"_nucl.db -parse_seqids -out "${assembly%.*}"_nucl.db
done

#Search the baits sequences and NLR genes against the databases of sequences of Ren-Seq reads/assemblies.
for db in barcode11_emily_trimmed_all.trimmedReads_nucl.db barcode12_fenella_trimmed_all.trimmedReads_nucl.db lorma_barcode11_emily_all_nucl.db lorma_barcode12_fenella_all_nucl.db smartdenovo_barcode11_emily_trimmed_all_racon_round_10_nucl.db smartdenovo_barcode12_fenella_trimmed_all_racon_round_10_nucl.db
do
blastn -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand"  -num_threads 1 -max_target_seqs 1000000000 -evalue 0.0000000001 -query probes-R4-final.fas -db $db >> probes-R4-final.fas_vs_$db
qsub $scripts/sub_blastn_renseq.sh vesca_v1.1_nblrrs_augustus_cds.fasta $db
done

#Generate a file with sequence lengths for CDS input files.
#Note for future reference: BLAST can automatically output these values if fields slen and qlen specified.
for a in *.fasta
do
python $scripts/write_seq_length.py $a
done

for sample in barcode11_emily_trimmed_all.trimmedReads barcode12_fenella_trimmed_all.trimmedReads lorma_barcode11_emily_all lorma_barcode12_fenella_all smartdenovo_barcode11_emily_trimmed_all_racon_round_10 smartdenovo_barcode12_fenella_trimmed_all_racon_round_10
do
#Identify RBB matches.
python /home/sobczm/bin/popgen/other/rbb.py ${sample}.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db vesca_v1.1_nblrrs_augustus_cds.fasta_vs_${sample}_nucl.db >${sample}_vesca_v1.1_nblrrs_augustus_cds_rbb
#Count up the length of CDS, 5' UTR, 3' UTR using the original renseq vs vesca CDS BLAST output
python $scripts/blast_matches_summary2.py vesca_v1.1_nblrrs_augustus_cds_lengths.txt ${sample}_lengths.txt ${sample}.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db > ${sample}.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_cds
#Filter those results to only retain the RBB hits.
python $scripts/filter_by_two_columns.py ${sample}_vesca_v1.1_nblrrs_augustus_cds_rbb ${sample}.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_cds > ${sample}.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_cds_rbb_filt
python $scripts/filter_by_two_columns.py ${sample}_vesca_v1.1_nblrrs_augustus_cds_rbb ${sample}.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_5utr > ${sample}.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_5utr_rbb_filt 
python $scripts/filter_by_two_columns.py ${sample}_vesca_v1.1_nblrrs_augustus_cds_rbb ${sample}.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_3utr > ${sample}.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_3utr_rbb_filt 
done

#Plot figures summarising the results above in R.
Rscript --vanilla $scripts/stats_cds.R barcode11_emily_trimmed_all.trimmedReads.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_cds_rbb_filt barcode12_fenella_trimmed_all.trimmedReads.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_cds_rbb_filt EMFE_canu_vs_vesca_v1.1_nblrrs_augustus_cds

Rscript --vanilla $scripts/stats_cds.R smartdenovo_barcode11_emily_trimmed_all_racon_round_10.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_cds_rbb_filt smartdenovo_barcode12_fenella_trimmed_all_racon_round_10.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_cds_rbb_filt EMFE_smartdenovo_polished_assembly.fasta_vs_vesca_v1.1_nblrrs_augustus_cds

Rscript --vanilla $scripts/stats_cds.R lorma_barcode11_emily_all.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_cds_rbb_filt lorma_barcode12_fenella_all.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_cds_rbb_filt EMFE_lorma_polished_assembly.fasta_vs_vesca_v1.1_nblrrs_augustus_cds

for utr in 3 5
do
Rscript --vanilla $scripts/stats_utr.R barcode11_emily_trimmed_all.trimmedReads.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_${utr}utr_rbb_filt barcode12_fenella_trimmed_all.trimmedReads.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_${utr}utr_rbb_filt EMFE_canu_vs_vesca_v1.1_nblrrs_augustus_$utr $utr
Rscript --vanilla $scripts/stats_utr.R smartdenovo_barcode11_emily_trimmed_all_racon_round_10.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_${utr}utr_rbb_filt smartdenovo_barcode12_fenella_trimmed_all_racon_round_10.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_${utr}utr_rbb_filt EMFE_smartdenovo_vs_vesca_v1.1_nblrrs_augustus_$utr $utr
Rscript --vanilla $scripts/stats_utr.R lorma_barcode11_emily_all.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_${utr}utr_rbb_filt lorma_barcode12_fenella_all.fasta_vs_vesca_v1.1_nblrrs_augustus_cds_nucl.db_${utr}utr_rbb_filt EMFE_lorma_vs_vesca_v1.1_nblrrs_augustus_$utr $utr
done

#Convert all to png
for my_pdf in *.pdf
do
convert -verbose -density 500 "${my_pdf}" "${my_pdf%.*}.png"
done

#Quick analysis of the NLR parser output.
cd $input/analysis/NLRParser
for my_table in *.tsv
do
echo $my_table >>nlr_summary.txt
echo "complete" >>nlr_summary.txt
cat $my_table | awk '$3=="complete" {print $0}' | wc -l >>nlr_summary.txt
echo "partial" >>nlr_summary.txt
cat $my_table | awk '$3=="partial" {print $0}' | wc -l >>nlr_summary.txt
echo "pseudogene" >>nlr_summary.txt
cat $my_table | awk '$3=="pseudogene" {print $0}' | wc -l >>nlr_summary.txt
done
