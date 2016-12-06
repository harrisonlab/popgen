#!/bin/bash
#Converts from VCF to the structure file format.
#Takes two arguments:
#First: haploid|diploid
#Second: input file
#Example usage: sh prepare_structure.sh haploid test.vcf

ploidy=$1
input=$2

filename=$(basename "$input")
outfile="${filename%.*}.struc"

pgdspid=/home/sobczm/bin/PGDSpider_2.1.0.3

cd $PWD

if [ "$ploidy" = "haploid" ]
then
    conver="vcf_to_structure_haploid.spid"
elif [ "$ploidy" = "diploid" ]
then
    conver="vcf_to_structure_diploid.spid"
else
    echo "Error, wrong first argument"
fi

java -jar $pgdspid/PGDSpider2-cli.jar -inputfile $input \
-inputformat VCF -outputfile $outfile -outputformat STRUCTURE -spid $pgdspid/$conver

dos2unix $outfile
