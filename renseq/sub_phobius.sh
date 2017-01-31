#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l virtual_free=1G
#$ -l h=blacklace02.blacklace|blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace


#-----
# 1. Set variables
#-----

CurPath=$PWD
InFile="$1"

Proteins=$(basename $InFile)
echo "$InFile"

phobius.pl $Infile > ${Infile%.fa*}_phoebius.txt
