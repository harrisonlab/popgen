#!/bin/bash
scripts=/home/sobczm/bin/popgen/snp
input=/home/sobczm/popgen/snp/snp_chip/crown_rot_gwas
cd $input

#First, running plink and TASSEL without imputation of genotypes - as we don't really know if we can impute them correctly, anyway - same problem as phasing.