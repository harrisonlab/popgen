
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