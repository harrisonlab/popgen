#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 2
#$ -l virtual_free=1G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace11.blacklace|blacklace12.blacklace

#Add '--prior logistic' when invoking structure below to change prior for more sensitive but exponentially slower structure detection

input=$1
k=$2

if [ -z "$3" ]
then
    model=simple
else
    model=logistic
fi

structure=/home/sobczm/bin/fastStructure
python $structure/structure.py -K $k --input $input --output $input --prior $model
