#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l virtual_free=1G
#$ -l h=blacklace02.blacklace|blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

#Cluster Ren-Seq baits to avoid redundancy. Input arguments:
#First - input fasta file with bait design
#Second - Max fraction identity allowed

fasta=$1
id_threshold=$2

filename=$(basename "$fasta")
name=${filename%.*}

usearch=/home/sobczm/bin/usearch/usearch9.0.2132_i86linux32

$usearch -cluster_fast $fasta -id $id_threshold -sort length -consout $"{name}"_$"{id_threshold}"_clust.fasta -uc $"{name}"_$"{id_threshold}".clusters
