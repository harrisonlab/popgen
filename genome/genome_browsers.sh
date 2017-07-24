#!/bin/bash
#This carried out on the genome server.

#Install UCSC and Ensembl genome browsers for ananassa latest genome.
#USCS tools (kentUtils) installed as well GBiB Assembly Hubs.
#Aim: create an assembly hub for the latest version of the ananassa genome along with
#its annotation by D. Swoerbreck group.
input=/home/sobczm/ananassa_ah
cp /home/sobczm/ananassa_annotation/Data_for_Maria_Sobczyk_02May2017/genome/v1_ns_remapped.fasta $input
#The .2bit file is constructed from the fasta sequence for the assembly. 
faToTwoBit v1_ns_remapped.fasta ananassa_05_17.2bit
#Use the twoBitInfo to verify the sequences in this assembly and create a chrom.sizes file which is not used in the hub, but is useful in later processing to construct the big* files:
twoBitInfo ananassa_05_17.2bit stdout | sort -k2rn > ananassa_05_17.chrom.sizes