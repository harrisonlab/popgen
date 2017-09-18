#!/bin/bash
scripts=/home/sobczm/bin/popgen/renseq
input=/home/sobczm/popgen/renseq/strawberry/reads

#Copy raw data over to start the analysis.
cp -r /data/seq_data/minion/2017/20170823_1707_RENseq-Em-Fen/20170823_1714_RENseq-Em-Fen/albacore_output_1.2.4 $input

#Extract the reads with porechop


#Blast the raw reads to establish % hits R genes.