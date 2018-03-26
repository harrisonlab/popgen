#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 6
#$ -l virtual_free=1G
#$ -l h=blacklace02.blacklace|blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

source /home/armita/.profile

DirIn=$1
FastqOut=$2

#nanopolish index -d $DirIn $FastqOut
CurrDir=$PWD
cd $DirIn
for a in $(echo ls -ld */)
do
nanopolish extract --type template $DirIn/$a >>$CurrDir/$FastqOut
done