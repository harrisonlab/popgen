#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l virtual_free=1G
#$ -l h=blacklace02.blacklace|blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

#set -u
#set -e

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

WorkDir=$TMPDIR/blastall
mkdir -p $WorkDir
cd $WorkDir

for File in $(ls "$BlastDB"*); do
  cp $File .
done
cp $Infile $Query

echo "hits will be named:"
echo $Hits
echo "hits will be moved to:"
echo "$OutFile"

blastall -d $Db -p blastp -i $Query -v 100000 -b 100000 -e 1e-5 -m 8 -F 'm S' -a 1 -o $Hits
cp $Hits $OutFile
