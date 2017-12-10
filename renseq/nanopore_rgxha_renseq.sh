
##Re-call 1D^2 reads using Albacore on the nanopore node.
~/.local/bin/full_1dsq_basecaller.py \
--flowcell FLO-MIN107 \
--kit SQK-LSK308 \
--input /home/nanopore/20171030_gridion_Redgauntlet_Hapil-RENseq/reads \
--recursive \
--worker_threads 12 \
--save_path /home/nanopore/20171030_gridion_Redgauntlet_Hapil-RENseq/reads/recalled \
--output_format fastq,fast5 \
--reads_per_fastq_batch 4000

scripts=/home/sobczm/bin/popgen/renseq
input=/home/sobczm/popgen/renseq/strawberry/reads/albacore_rgxha
cd $input
#Copy the pass reads across from the 1d^2 protocol. As few reads are "pass", going to recall using the 1d protocol.
scp -r nanopore@nanopore:/home/nanopore/20171030_gridion_Redgauntlet_Hapil-RENseq/reads/recalled/1dsq_analysis/workspace/pass ./

~/.local/bin/read_fast5_basecaller.py \
--flowcell FLO-MIN107 \
--kit SQK-LSK108 \
--input /data/scratch/nanopore_tmp_data/sobczm/reads \
--recursive \
--worker_threads 12 \
--save_path /home/nanopore/20171030_gridion_Redgauntlet_Hapil-RENseq/1d \
--output_format fastq,fast5 \
--reads_per_fastq_batch 4000

#Trim the adapters. Separate out RG and HA samples.
input=/home/sobczm/popgen/renseq/strawberry/reads/albacore_rgxha
scripts=/home/sobczm/bin/popgen/renseq
cd $input

#In this analysis, will combine 1d2 'pass' reads with 1d 'pass' reads, and in the case of reads called as 'pass' in those, will prioritise 1d2 over 1d.

#Concatenate reads: 1D2
for f in $input/pass/*.fastq
do
cat $f >>$input/rgxha_1d2_reads.fastq
done

#Concatenate reads: 1D
for f in $input/20171030_gridion_Redgauntlet_Hapil-RENseq/1d/workspace/pass/*.fastq
do
cat $f >>$input/rgxha_1d_reads.fastq
done

python $scripts/combine_1d_1d2_pass.py rgxha_1d2_reads.fastq rgxha_1d_reads.fastq >rgxha_1d2_1d_combined.fastq

#Run unmodified porechop on those reads
qsub $scripts/sub_poretools_demultiplex.sh rgxha_1d2_1d_combined.fastq.gz porechop_native
#As expected, no Nanopore barcodes detected.
#Now, running modified porechop with custom Illumina barcodes using modified porechop version from https://github.com/Psy-Fer/Porechop
#As the Illumina barcodes differ only by 5 nucleotides it is too difficult for reads to be binned correctly - Nanopore barcodes are designed to maximize sequence differences.
#The modified adapter file used for that is located in /home/sobczm/bin/porechop_alt/Porechop/alt_adapters.py