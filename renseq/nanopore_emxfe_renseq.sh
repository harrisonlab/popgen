#!/bin/bash
scripts=/home/sobczm/bin/popgen/renseq
input=/home/sobczm/popgen/renseq/strawberry/reads/albacore_emxfe
cd $input
#Concatenate Extract the confirmed reads with porechop (Barcode 11 and Barcode 12 only now)
for f in /data/seq_data/minion/2017/20170823_1707_RENseq-Em-Fen/20170823_1714_RENseq-Em-Fen/albacore_output_1.2.4/workspace/barcode11/*.fastq
do
cat $f >>$input/barcode11_emily.fastq
done

for f in /data/seq_data/minion/2017/20170823_1707_RENseq-Em-Fen/20170823_1714_RENseq-Em-Fen/albacore_output_1.2.4/workspace/barcode12/*.fastq
do
cat $f >>$input/barcode12_fenella.fastq
done


#Convert to FASTA
/home/sobczm/bin/seqtk/seqtk seq -a barcode11_emily.fastq > barcode11_emily.fasta
/home/sobczm/bin/seqtk/seqtk seq -a barcode12_fenella.fastq > barcode12_fenella.fasta


#What is read length?
perl /home/sobczm/bin/popgen/other/assemblathon_stats.pl barcode11_emily.fasta
     Number of scaffolds     627300
                                     Total size of scaffolds 3321273174
                                            Longest scaffold      50608
                                           Shortest scaffold        128
                                 Number of scaffolds > 1K nt     620576  98.9%
                                Number of scaffolds > 10K nt       2080   0.3%
                               Number of scaffolds > 100K nt          0   0.0%
                                 Number of scaffolds > 1M nt          0   0.0%
                                Number of scaffolds > 10M nt          0   0.0%
                                          Mean scaffold size       5295
                                        Median scaffold size       5291
                                         N50 scaffold length       5366
                                          L50 scaffold count     276441
perl /home/sobczm/bin/popgen/other/assemblathon_stats.pl barcode12_fenella.fasta
          Number of scaffolds     600337
                                     Total size of scaffolds 3082283572
                                            Longest scaffold      22473
                                           Shortest scaffold        158
                                 Number of scaffolds > 1K nt     591240  98.5%
                                Number of scaffolds > 10K nt       2053   0.3%
                               Number of scaffolds > 100K nt          0   0.0%
                                 Number of scaffolds > 1M nt          0   0.0%
                                Number of scaffolds > 10M nt          0   0.0%
                                          Mean scaffold size       5134
                                        Median scaffold size       5165
                                         N50 scaffold length       5267
                                          L50 scaffold count     259702
                                                 scaffold %A      29.81
                                                 scaffold %C      19.76
                                                 scaffold %G      20.02
                                                 scaffold %T      30.41
                                                 scaffold %N       0.00
                                         scaffold %non-ACGTN       0.00
                             Number of scaffold non-ACGTN nt          0

sh /home/sobczm/bin/popgen/other/run_fastqc.sh barcode11_emily.fastq
cat barcode11_emily.fastq | awk '{if(NR%4==2) print length($1)}' | sort -n > emily_read_length.txt
sh /home/sobczm/bin/popgen/other/run_fastqc.sh barcode12_fenella.fastq
cat barcode12_fenella.fastq | awk '{if(NR%4==2) print length($1)}' | sort -n  > fenella_read_length.txt

#Blast the raw reads to establish % hits R genes.
#Copy the blast db
cp -r /home/sobczm/popgen/renseq/strawberry/reads/Helen_Bates_EMR.RH.ENQ-1704.A.01/analysis/vesca_v1.1_nblrrs_augustus_mrna_nucl.db* ./

for ass in $input/barcode11_emily.fasta $input/barcode12_fenella.fasta
do
    blastn -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand"  -num_threads 1 -max_target_seqs 1 -evalue 0.0000000001 -query $ass -db vesca_v1.1_nblrrs_augustus_mrna_nucl.db >> $(basename $ass)_vs_$db
done

#How many reads have matches?
 cut -f1 barcode11_emily.fasta_vs_vesca_v1.1_nblrrs_augustus_mrna_nucl.db | sort | uniq | wc
 498443  498443 18442391

 cut -f1 barcode12_fenella.fasta_vs_vesca_v1.1_nblrrs_augustus_mrna_nucl.db | sort | uniq | wc
 478224  478224 17694288

#How many have matches at least 1 kbp long?
awk -F"\t" '$4 > 999 { print $1 }' barcode11_emily.fasta_vs_vesca_v1.1_nblrrs_augustus_mrna_nucl.db | sort | uniq | wc
418522  418522 15485314
awk -F"\t" '$4 > 999 { print $1 }' barcode12_fenella.fasta_vs_vesca_v1.1_nblrrs_augustus_mrna_nucl.db | sort | uniq | wc
 395610  395610 14637570

#Demultiplex the unclassified reads (permissive search settings) to see if any barcode11 or barcode12 sequences found there.
#Concatenate and gzip all FASTQ files in the unclassified bin.
porechop=/home/armita/prog/porechop/Porechop
albacoredir=/data/seq_data/minion/2017/20170823_1707_RENseq-Em-Fen/20170823_1714_RENseq-Em-Fen/albacore_output_1.2.4/workspace
for fast in $albacoredir/unclassified/*.fastq
do
cat $fast >>all_unclassified.fastq
done
gzip all_unclassified.fastq

#Demultiplex and trim adapters from unclassified reads using porechop
qsub $scripts/sub_poretools_demultiplex.sh all_unclassified.fastq.gz porechop_unclassified
#Retrieve additional barcode 11 reads
gzip -d $input/porechop_unclassified/NB11.fastq.gz
#Retrive additional barcode 12 reads
gzip -d $input/porechop_unclassified/NB12.fastq.gz

#Convert to FASTA
/home/sobczm/bin/seqtk/seqtk seq -a $input/porechop_unclassified/NB11.fastq > $input/porechop_unclassified/NB11.fasta
/home/sobczm/bin/seqtk/seqtk seq -a $input/porechop_unclassified/NB12.fastq > $input/porechop_unclassified/NB12.fasta

#Split reads and trim adapters using porechop - barcode 11
gzip barcode11_emily.fastq
qsub $scripts/sub_poretools_trim.sh barcode11_emily.fastq.gz
gzip -d barcode11_emily_trim.fastq.gz 
/home/sobczm/bin/seqtk/seqtk seq -a barcode11_emily_trim.fastq > barcode11_emily_trim.fasta

#Split reads and trim adapters using porechop - barcode 12
gzip barcode12_fenella.fastq
qsub $scripts/sub_poretools_trim.sh barcode12_fenella.fastq.gz
gzip -d barcode12_fenella_trim.fastq.gz
/home/sobczm/bin/seqtk/seqtk seq -a barcode12_fenella_trim.fastq > barcode12_fenella_trim.fasta
#Two rounds of nanocorrect of trimmed reads.
#First concatenate additional (previously unclassified) reads with the main trimmed reads
cat barcode11_emily_trim.fasta >> barcode11_emily_trimmed_all.fasta
cat $input/porechop_unclassified/NB11.fasta >> barcode11_emily_trimmed_all.fasta

cat barcode12_fenella_trim.fasta >> barcode12_fenella_trimmed_all.fasta
cat $input/porechop_unclassified/NB12.fasta >> barcode12_fenella_trimmed_all.fasta

cat barcode11_emily_trim.fastq >> barcode11_emily_trimmed_all.fastq
cat $input/porechop_unclassified/NB11.fastq >> barcode11_emily_trimmed_all.fastq

cat barcode12_fenella_trim.fastq >> barcode12_fenella_trimmed_all.fastq
cat $input/porechop_unclassified/NB12.fastq >> barcode12_fenella_trimmed_all.fastq


###Generate a dataset with random 75% of reads and put them through all the steps below to check the effect of a lower number of reads on the R gene detection.
/home/sobczm/bin/seqtk/seqtk sample barcode11_emily_trimmed_all.fastq 450000 > barcode11_emily_trimmed_all_075.fastq
/home/sobczm/bin/seqtk/seqtk sample barcode12_fenella_trimmed_all.fastq 450000 > barcode12_fenella_trimmed_all_075.fastq
/home/sobczm/bin/seqtk/seqtk sample barcode11_emily_trimmed_all.fasta 450000 > barcode11_emily_trimmed_all_075.fasta
/home/sobczm/bin/seqtk/seqtk sample barcode12_fenella_trimmed_all.fasta 450000 > barcode12_fenella_trimmed_all_075.fasta
#Generate a polished set of reads - two rounds of error correcting should lead to 97% accuracy.
qsub $scripts/sub_nanocorrect.sh barcode11_emily_trimmed_all.fasta nanocorrect_emily_nanocorrect
qsub $scripts/sub_nanocorrect.sh barcode12_fenella_trimmed_all.fasta nanocorrect_fenella_nanocorrect
qsub $scripts/sub_nanocorrect.sh barcode11_emily_trimmed_all_075.fasta nanocorrect_emily_nanocorrect_075
qsub $scripts/sub_nanocorrect.sh barcode12_fenella_trimmed_all_075.fasta nanocorrect_fenella_nanocorrect_075

#Cannot make it work, use LoRMA instead.
qsub $scripts/sub_lorma.sh barcode11_emily_trimmed_all.fasta  lorma_barcode11_emily_all
qsub $scripts/sub_lorma.sh barcode12_fenella_trimmed_all.fasta lorma_barcode12_fenella_all
qsub $scripts/sub_lorma.sh barcode11_emily_trimmed_all_075.fasta lorma_barcode11_emily_075
qsub $scripts/sub_lorma.sh barcode12_fenella_trimmed_all_075.fasta lorma_barcode12_fenella_075

#Read error correction with Canu
gzip barcode11_emily_trimmed_all.fastq 
qsub $scripts/sub_canu_correct.sh barcode11_emily_trimmed_all.fastq.gz 25m barcode11_emily_trimmed_all barcode11_emily_trimmed_all

gzip barcode12_fenella_trimmed_all.fastq 
gzip barcode11_emily_trimmed_all_075.fastq 
gzip barcode12_fenella_trimmed_all_075.fastq 
qsub $scripts/sub_canu_correct.sh barcode12_fenella_trimmed_all.fastq.gz 25m barcode12_fenella_trimmed_all barcode12_fenella_trimmed_all
qsub $scripts/sub_canu_correct.sh barcode11_emily_trimmed_all_075.fastq.gz 25m barcode11_emily_trimmed_all barcode11_emily_trimmed_all_075
qsub $scripts/sub_canu_correct.sh barcode12_fenella_trimmed_all_075.fastq.gz 25m barcode12_fenella_trimmed_all barcode12_fenella_trimmed_all_075

#Assembly of FASTQ reads with smartdenovo
qsub $scripts/sub_SMARTdenovo.sh barcode12_fenella_trimmed_all/barcode12_fenella_trimmed_all.trimmedReads.fasta.gz smartdenovo_barcode12_fenella_trimmed_all smartdenovo_barcode12_fenella_trimmed_all
qsub $scripts/sub_SMARTdenovo.sh barcode12_fenella_trimmed_all_075/barcode12_fenella_trimmed_all.trimmedReads.fasta.gz smartdenovo_barcode12_fenella_trimmed_075 smartdenovo_barcode12_fenella_trimmed_075
qsub $scripts/sub_SMARTdenovo.sh barcode11_emily_trimmed_all/barcode11_emily_trimmed_all.trimmedReads.fasta.gz smartdenovo_barcode11_emily_trimmed_all smartdenovo_barcode11_emily_trimmed_all
qsub $scripts/sub_SMARTdenovo.sh barcode11_emily_trimmed_all_075/barcode11_emily_trimmed_all.trimmedReads.fasta.gz smartdenovo_barcode11_emily_trimmed_075 smartdenovo_barcode11_emily_trimmed_075
#Error correction using racon
ass=smartdenovo_barcode12_fenella_trimmed_all 
qsub $scripts/sub_racon.sh $ass/$ass.dmo.lay.utg barcode12_fenella_trimmed_all.fastq.gz 10 $ass/racon
ass=smartdenovo_barcode12_fenella_trimmed_075 
qsub $scripts/sub_racon.sh $ass/$ass.dmo.lay.utg barcode12_fenella_trimmed_all_075.fastq.gz 10 $ass/racon
ass=smartdenovo_barcode11_emily_trimmed_all 
qsub $scripts/sub_racon.sh $ass/$ass.dmo.lay.utg barcode11_emily_trimmed_all.fastq.gz 10 $ass/racon
ass=smartdenovo_barcode11_emily_trimmed_075
qsub $scripts/sub_racon.sh $ass/$ass.dmo.lay.utg barcode11_emily_trimmed_all.fastq.gz 10 $ass/racon

#Assembly correction using nanopolish
#Re-extract reads
qsub $scripts/sub_nanopolish_extract.sh $raw_reads/barcode11 nanopolish_barcode11_emily.fasta
qsub $scripts/sub_nanopolish_extract.sh $raw_reads/barcode12 nanopolish_barcode12_fenella.fasta
gzip nanopolish_barcode11_emily.fasta
gzip nanopolish_barcode12_fenella.fasta 
#Align reads
qsub $scripts/sub_bwa_nanopolish.sh smartdenovo_barcode11_emily_trimmed_075/racon/smartdenovo_barcode11_emily_trimmed_075_racon_round_10.fasta nanopolish_barcode11_emily.fasta.gz smartdenovo_barcode11_emily_trimmed_075/nanopolish
qsub $scripts/sub_bwa_nanopolish.sh smartdenovo_barcode11_emily_trimmed_all/racon/smartdenovo_barcode11_emily_trimmed_all_racon_round_10.fasta nanopolish_barcode11_emily.fasta.gz smartdenovo_barcode11_emily_trimmed_all/nanopolish
qsub $scripts/sub_bwa_nanopolish.sh smartdenovo_barcode12_fenella_trimmed_075/racon/smartdenovo_barcode12_fenella_trimmed_075_racon_round_10.fasta nanopolish_barcode12_fenella.fasta.gz smartdenovo_barcode12_fenella_trimmed_075/nanopolish
qsub $scripts/sub_bwa_nanopolish.sh smartdenovo_barcode12_fenella_trimmed_all/racon/smartdenovo_barcode12_fenella_trimmed_all_racon_round_10.fasta nanopolish_barcode12_fenella.fasta.gz smartdenovo_barcode12_fenella_trimmed_all/nanopolish
#Submit alignments for nanopolish - variant calling.
NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
python $NanoPolishDir/nanopolish_makerange.py smartdenovo_barcode11_emily_trimmed_075/racon/smartdenovo_barcode11_emily_trimmed_075_racon_round_10.fasta > smartdenovo_barcode11_emily_trimmed_075/nanopolish/nanopolish_range.txt
python $NanoPolishDir/nanopolish_makerange.py smartdenovo_barcode11_emily_trimmed_all/racon/smartdenovo_barcode11_emily_trimmed_all_racon_round_10.fasta > smartdenovo_barcode11_emily_trimmed_all/nanopolish/nanopolish_range.txt
python $NanoPolishDir/nanopolish_makerange.py smartdenovo_barcode12_fenella_trimmed_075/racon/smartdenovo_barcode12_fenella_trimmed_075_racon_round_10.fasta > smartdenovo_barcode12_fenella_trimmed_075/nanopolish/nanopolish_range.txt
python $NanoPolishDir/nanopolish_makerange.py smartdenovo_barcode12_fenella_trimmed_all/racon/smartdenovo_barcode12_fenella_trimmed_all_racon_round_10.fasta > smartdenovo_barcode12_fenella_trimmed_all/nanopolish/nanopolish_range.txt

Assembly=smartdenovo_barcode11_emily_trimmed_075/racon/smartdenovo_barcode11_emily_trimmed_075_racon_round_10.fasta 
RawReads=nanopolish_barcode11_emily.fasta.gz 
AlignedReads=smartdenovo_barcode11_emily_trimmed_075/nanopolish/reads.sorted.bam
Ploidy=2
OutDir=smartdenovo_barcode11_emily_trimmed_075/nanopolish
for Region in $(cat $OutDir/nanopolish_range.txt | tail -n+21); do
Jobs=$(qstat | grep 'sub_nanopo' | wc -l)
while [ $Jobs -gt 4 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_nanopo' | wc -l)
done		
printf "\n"
echo $Region
echo $Region >> nanopolish_log.txt
qsub $scripts/sub_nanopolish_variants.sh $Assembly $RawReads $AlignedReads $Ploidy $Region $OutDir/$Region
done

for variants in $OutDir/*/*.txt
do
cat $variants >> $OutDir/nanopolish_all_variants.txt
done

Assembly=smartdenovo_barcode11_emily_trimmed_all/racon/smartdenovo_barcode11_emily_trimmed_all_racon_round_10.fasta 
RawReads=nanopolish_barcode11_emily.fasta.gz 
AlignedReads=smartdenovo_barcode11_emily_trimmed_all/nanopolish/reads.sorted.bam
Ploidy=2
OutDir=smartdenovo_barcode11_emily_trimmed_all/nanopolish

for Region in $(cat $OutDir/nanopolish_range.txt | tail -n+21); do
Jobs=$(qstat | grep 'sub_nanopo' | wc -l)
while [ $Jobs -gt 4 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_nanopo' | wc -l)
done		
printf "\n"
echo $Region
echo $Region >> nanopolish_log.txt
qsub $scripts/sub_nanopolish_variants.sh $Assembly $RawReads $AlignedReads $Ploidy $Region $OutDir/$Region
done

for variants in $OutDir/*/*.txt
do
cat $variants >> $OutDir/nanopolish_all_variants.txt
done

Assembly=smartdenovo_barcode12_fenella_trimmed_075/racon/smartdenovo_barcode12_fenella_trimmed_075_racon_round_10.fasta 
RawReads=nanopolish_barcode12_fenella.fasta.gz 
AlignedReads=smartdenovo_barcode12_fenella_trimmed_075/nanopolish/reads.sorted.bam
Ploidy=2
OutDir=smartdenovo_barcode12_fenella_trimmed_075/nanopolish

for Region in $(cat $OutDir/nanopolish_range.txt | tail -n+21); do
Jobs=$(qstat | grep 'sub_nanopo' | wc -l)
while [ $Jobs -gt 4 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_nanopo' | wc -l)
done		
printf "\n"
echo $Region
echo $Region >> nanopolish_log.txt
qsub $scripts/sub_nanopolish_variants.sh $Assembly $RawReads $AlignedReads $Ploidy $Region $OutDir/$Region
done

for variants in $OutDir/*/*.txt
do
cat $variants >> $OutDir/nanopolish_all_variants.txt
done

Assembly=smartdenovo_barcode12_fenella_trimmed_all/racon/smartdenovo_barcode12_fenella_trimmed_all_racon_round_10.fasta 
RawReads=nanopolish_barcode12_fenella.fasta.gz 
AlignedReads=smartdenovo_barcode12_fenella_trimmed_all/nanopolish/reads.sorted.bam
Ploidy=2
OutDir=smartdenovo_barcode12_fenella_trimmed_all/nanopolish

for Region in $(cat $OutDir/nanopolish_range.txt | tail -n+21); do
Jobs=$(qstat | grep 'sub_nanopo' | wc -l)
while [ $Jobs -gt 4 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_nanopo' | wc -l)
done		
printf "\n"
echo $Region
echo $Region >> nanopolish_log.txt
qsub $scripts/sub_nanopolish_variants.sh $Assembly $RawReads $AlignedReads $Ploidy $Region $OutDir/$Region
done

for variants in $OutDir/*/*.txt
do
cat $variants >> $OutDir/nanopolish_all_variants.txt
done

#Align Illumina reads and carry out variant calling with GATK.