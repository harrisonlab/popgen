#!/bin/bash
input=/home/sobczm/popgen/snp/snp_chip/haplotypes/comparison
scripts=/home/sobczm/bin/popgen/snp

#Copy tha map-based phasing results.

#Extract results from matching individuals - haplotyping.
for a in *_map.ped
do
python $scripts/overall_haplotype_similarity.py $a ${a%_map.ped}_impute.ped ${a%.ped}.info ${a%_map.ped}_impute.info
python $scripts/short_identical_haplotypes.py $a ${a%_map.ped}_impute.ped ${a%.ped}.info ${a%_map.ped}_impute.info
python $scripts/haplotype_length_distribution.py $a ${a%_map.ped}_impute.ped ${a%.ped}.info ${a%_map.ped}_impute.info
done

#Plot results using ggplot2. 
for a in *_overall.stat
do
Rscript --vanilla $scripts/overall_haplotype_similarity_fig.R $a 8 4
done
