#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l virtual_free=1G
#$ -l h=blacklace02.blacklace|blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace


CurPath=$PWD
BlastDB=$1
Infile=$2
OutFile=$3

echo "Infiles:"
echo $BlastDB
echo $Infile
echo $OutFile

Db=$(basename $1)
Query=$(basename $2)
Hits=$(basename $3)

echo "Copying to Files:"
echo $Db
echo $Query

echo "Outfiles will be named:"
echo $Hits


echo "hits will be named:"
echo $Hits
echo "hits will be moved to:"
echo "$OutFile"

blastall -d $BlastDB -p blastp -i $Infile -v 100000 -b 100000 -e 1e-5 -m 8 -F 'm S' -a 1 -o $OutFile

echo $Db
echo $Query
echo $Hits
