#!/bin/bash
#Converts from FASTA to VCF (outputs only variable positions - SNPs)
#Takes one argumens:
#input file
#Example usage: sh fasta_to_vcf.sh test.vcf

input=$1

filename=$(basename "$input")
outfile="${filename%.*}.vcf"

pgdspid=/home/sobczm/bin/PGDSpider_2.1.0.3

cd $PWD

conver="fasta_to_vcf_snp.spid"

java -jar $pgdspid/PGDSpider2-cli.jar -inputfile $input \
-inputformat FASTA -outputfile $outfile -outputformat VCF -spid $pgdspid/$conver

dos2unix $outfile
