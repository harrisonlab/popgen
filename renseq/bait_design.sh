#!/bin/bash
input=/home/sobczm/popgen/renseq/input/transcriptomes
scripts=/home/sobczm/bin/popgen/renseq
usearch=/home/sobczm/bin/usearch/usearch9.0.2132_i86linux32

mkdir -p $input/baits
cp all_rgenes.fasta $input/baits
#Design baits based on selected R genes
n=4
python $scripts/create_baits.py --inp $input/baits/all_rgenes.fasta --coverage $n \
--out all_rgenes_baits_"$n"x.fasta
#Results in around 18100 baits. Needs to increase tilling density
n=5
#Results in around 22600 baits.
#Check for presence of Ns and repetitive sequences
fasta=all_rgenes_baits_"$n"x.fasta
#Hard-mask repetitive regions with Ns
$usearch -fastx_mask $fasta -qmask dust -fastaout "${fasta%.*}"_masked.fasta -hardmask
#Remove the baits containing Ns.
python $scripts/remove_N_fasta.py all_rgenes_baits_"$n"x_masked.fasta
#Results in retention of around 21775 baits, ie. around 4% baits discarded
#Cluster the baits at 95% identity
fasta=all_rgenes_baits_"$n"x_masked_noN.fasta
id_threshold=0.95
$usearch -cluster_fast $fasta -id $id_threshold -sort length \
-consout ${fasta%.*}_"$id_threshold"_clust.fasta -uc ${fasta%.*}_"$id_threshold".clusters
#Out of 21776 seqs, only 11825 are unique!

#Need to increase tilling density - this time around 8x.
n=8
id_threshold=0.95
#to fix the ends for no end gaps, use -leftjust and -rightjust
$usearch -cluster_fast $fasta -id $id_threshold -sort length -leftjust -rightjust \
-consout ${fasta%.*}_"$id_threshold"_clust.fasta -uc ${fasta%.*}_"$id_threshold".clusters
#... And get 18250 uniques

#Finally, trying out 12x tilling
n=12
id_threshold=0.95
#... And get 26129 uniques, as before clustering into 25067 clusters at 90% identity
# Clustering at 99% percent identity:
id_threshold=0.99
#Get similar results

#Lastly, do very high tilling coverage
n=20
id_threshold=0.99
#Get 40074 uniques

#Select 25,000 random sequences to be submitted to MYbaits for testing
