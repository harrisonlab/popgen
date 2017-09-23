#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 6
#$ -l virtual_free=2G
#$ -l h=blacklace02.blacklace|blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

lorma=/home/sobczm/bin/LoRMA-0.4
input_fasta=$1

$lorma/lorma.sh $input_fasta