#!/bin/bash
scripts=/home/sobczm/bin/popgen/renseq
input=/home/sobczm/popgen/renseq/apple

#Download all the apple genome resources
cd $input
mkdir genome
cd genome
wget -r ftp://climb.genomics.cn/pub/10.5524/100001_101000/100189/
#File with CDS sequences
ls -l $input/genome/EVM.out.cds
#File with pep sequences
ls -l $input/genome/gene/EVM.out.pep